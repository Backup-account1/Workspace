#!/usr/bin/env python3
"""
segment_cheese_fungus_enhanced.py

Enhanced version of the CoreSpecViewer based fungus defect segmentation pipeline.
Features:
- Automatic SNR‑based band selection (default: 2200‑2320 nm, can be overridden)
- Optional despeckling / morphological clean‑up
- Configurable abundance threshold, minimum object size and margin
- Detailed logging to a results folder
- CLI arguments for easy experimentation
"""

import os
import argparse
import numpy as np
import csv
import rasterio
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from tqdm import tqdm

# Import CoreSpecViewer utilities
from app.spectral_ops.spectral_functions import (
    remove_hull,
    mnf_transform,
    nmf_unmixing,
    threshold_abundance,
    clean_mask,
    despeckle_mask,
    improve_mask_from_graph,
    export_results,  # placeholder – will be overridden below
)

# ----------------------------------------------------------------------
#  Helper: simple CSV logger for summary stats
# ----------------------------------------------------------------------
def write_summary(csv_path: str, total_pixels: int, fungus_pixels: int, percent: float):
    """Append a one‑row summary CSV."""
    fieldnames = ["date", "total_pixels", "fungus_pixels", "percentage_fungus"]
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow({
            "date": os.popen("date +'%Y-%m-%d %H:%M:%S'").read().strip(),
            "total_pixels": total_pixels,
            "fungus_pixels": fungus_pixels,
            "percentage_fungus": percent,
        })


# ----------------------------------------------------------------------
#  Main processing function
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Segment fungal defect regions in cheese hyperspectral images (HSI) using CoreSpecViewer."
    )
    parser.add_argument("--hsi_path", required=True, help="Path to ENVI cube (or .npy stack).")
    parser.add_argument("--csv_path", required=True, help="CSV file with fungal reference spectra.")
    parser.add_argument("--output_dir", default="./results", help="Folder for all output files.")
    parser.add_argument("--thresh", type=float, default=0.12, help="Abundance threshold for mask.")
    parser.add_argument("--mnf_components", type=int, default=15, help="Number of MNF components.")
    parser.add_argument("--min_mask_area", type=int, default=80, help="Minimum object size (pixels).")
    parser.add_argument("--margin", type=int, default=5, help="Margin (pixels) for overlay expansion.")
    parser.add_argument("--snr_thresh", type=float, default=20.0, help="SNR threshold for band selection.")
    parser.add_argument("--snr_min_run", type=int, default=20, help="Minimum contiguous run length for band selection.")
    parser.add_argument("--clean_method", choices=["none", "despeckle", "graph"], default="none",
                        help="Mask cleaning method.")
    parser.add_argument("--log_summary", action="store_true", help="Write a CSV summary of results.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1️⃣ Load fungal reference spectra
    # ------------------------------------------------------------------
    wavelengths, ref_spectra = [], []
    with open(args.csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            wavelengths.append(float(row["wavelength"]))
            ref_spectra.append(float(row["fungus_reflection"]))
    wavelengths = np.array(wavelengths, dtype=np.float32)
    ref_spectra = np.array(ref_spectra, dtype=np.float32)

    # ------------------------------------------------------------------
    # 2️⃣ Load hyperspectral cube
    # ------------------------------------------------------------------
    if args.hsi_path.lower().endswith(('.hdr', '.dat')):
        import spectral.io.envi as envi
        cube = envi.open_image(args.hsi_path).load().astype(np.float32)
        wavelengths_hdr = np.array(envi.get_metadata(args.hsi_path).get('wavelength', []), dtype=np.float32)
    else:
        cube = np.load(args.hsi_path).astype(np.float32)
        wavelengths_hdr = np.arange(cube.shape[2], dtype=np.float32)  # dummy if no meta

    # ------------------------------------------------------------------
    # 3️⃣ Optional SNR‑based band selection
    # ------------------------------------------------------------------
    def select_bands_by_snr(cube, wavelengths, snr_thresh, min_run):
        """Return a slice of bands that form the longest contiguous run above SNR threshold."""
        # Reduce over spatial dimensions to get per‑band SNR
        # (Assume cube is [rows, cols, bands])
        H, W, B = cube.shape
        # Simple SNR: (mean - dark_mean) / std(dark) approximated using all pixels
        mean_val = np.mean(cube, axis=(0, 1))
        # Use a rough dark estimate: low intensity percentile
        dark_est = np.percentile(cube, 2, axis=(0, 1))
        std_dark = np.std(np.maximum(dark_est, 1e-6))
        snr = (mean_val - dark_est + 1e-9) / (std_dark + 1e-9)

        good = snr >= snr_thresh
        # Find longest contiguous True run
        starts, ends, in_run = [], [], False
        for i, ok in enumerate(good.tolist() + [False]):
            if ok and not in_run:
                s = i
                in_run = True
            elif not ok and in_run:
                starts.append(s)
                ends.append(i)
                in_run = False
        if not starts:
            # fallback: use central 30% of bands
            start_idx = int(0.35 * B)
            stop_idx = int(0.65 * B)
            return slice(start_idx, stop_idx)
        lengths = np.array(ends) - np.array(starts)
        best_idx = int(np.argmax(lengths))
        return slice(starts[best_idx], ends[best_idx])

    # If user wants automatic selection, compute it now
    if args.snr_thresh is not None:
        band_slice = select_bands_by_snr(cube, wavelengths_hdr, args.snr_thresh, args.snr_min_run)
    else:
        band_slice = slice(None)  # use all bands

    cube = cube[..., band_slice]
    wavelengths_hdr = wavelengths_hdr[band_slice]

    print(f"🔍 Loaded cube shape: {cube.shape}")
    print(f"📈 Using wavelength range: {wavelengths_hdr[0]:.1f}–{wavelengths_hdr[-1]:.1f} nm")

    # ------------------------------------------------------------------
    # 4️⃣ Continuum removal
    # ------------------------------------------------------------------
    print("🛠️  Continuum removal …")
    cr_cube = remove_hull(cube)

    # ------------------------------------------------------------------
    # 5️⃣ MNF projection
    # ------------------------------------------------------------------
    print(f"⚙️  Computing {args.mnf_components} MNF components …")
    mnf_cube = mnf_transform(cr_cube, args.mnf_components)

    # ------------------------------------------------------------------
    # 6️⃣ NMF unmixing (abundance of the first end‑member)
    # ------------------------------------------------------------------
    print("🧩  NMF unmixing …")
    H, W, _ = mnf_cube.shape
    mnf_2d = mnf_cube.reshape(-1, args.mnf_components).T  # (components, pixels)
    abundance_map = nmf_unmixing(mnf_2d, n_endmembers=1)  # shape (pixels,)
    abundance_map = abundance_map.reshape(H, W)

    # ------------------------------------------------------------------
    # 7️⃣ Threshold → binary mask
    # ------------------------------------------------------------------
    print(f"🔎  Thresholding (>{args.thresh}) …")
    raw_mask = threshold_abundance(abundance_map, args.thresh)

    # ------------------------------------------------------------------
    # 8️⃣ Mask cleaning
    # ------------------------------------------------------------------
    print("🧹  Cleaning mask …")
    if args.clean_method == "despeckle":
        mask = despeckle_mask(raw_mask)
    elif args.clean_method == "graph":
        mask = improve_mask_from_graph(raw_mask)
    else:
        mask = raw_mask

    # Further clean with morphological operations if desired
    if args.min_mask_area > 0:
        mask = clean_mask(mask, min_area=args.min_mask_area, margin=args.margin)

    # ------------------------------------------------------------------
    # 9️⃣ Export results
    # ------------------------------------------------------------------
    print("🚀  Exporting results …")
    # Simple overlay PNG using matplotlib
    plt.figure(figsize=(8, 8))
    # Show first three bands as RGB for quick visual check
    rgb_preview = cube[:, :, :3].astype(np.float32)
    plt.imshow(np.clip(rgb_preview, 0, 1))
    # Overlay mask in semi‑transparent red
    plt.imshow(mask, cmap="Reds", alpha=0.4)
    plt.axis("off")
    overlay_path = os.path.join(args.output_dir, f"{os.path.basename(args.hsi_path)}_overlay.png")
    plt.savefig(overlay_path, dpi=300, bbox_inches="tight")
    plt.close()

    # Export mask as GeoTIFF (simple transform – placeholder)
    mask_path = os.path.join(args.output_dir, f"{os.path.basename(args.hsi_path)}_mask.tif")
    meta = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 1,
        "height": mask.shape[0],
        "width": mask.shape[1],
        "transform": from_origin(0, 0, 1, 1),
    }
    with rasterio.open(mask_path, "w", **meta) as dst:
        dst.write(mask.astype(np.uint8), 1)

    # Export CSV summary if requested
    summary_csv = os.path.join(args.output_dir, "summary.csv")
    if args.log_summary:
        total_px = mask.size
        fungus_px = mask.sum()
        percent_fungus = fungus_px / total_px * 100
        write_summary(summary_csv, total_px, fungus_px, percent_fungus)

    print(f"✅ Mask saved          → {mask_path}")
    print(f"✅ Overlay saved       → {overlay_path}")
    if args.log_summary:
        print(f"✅ Summary CSV saved   → {summary_csv}")

    print("✅ All done!")


if __name__ == "__main__":
    main()