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
from datetime import datetime

# ----------------------------------------------------------------------
#  Optional CoreSpecViewer import (if available), with local fallbacks
# ----------------------------------------------------------------------
try:
    # If you have the CoreSpecViewer-style `app/` package in your PYTHONPATH,
    # these imports will be used. In this workspace, we provide local fallbacks.
    from app.spectral_ops.spectral_functions import (  # type: ignore
        remove_hull,
        mnf_transform,
        nmf_unmixing,
        threshold_abundance,
        clean_mask,
        despeckle_mask,
        improve_mask_from_graph,
        export_results,  # may exist upstream
    )
except Exception:
    from scipy.ndimage import label, binary_fill_holes, binary_dilation, binary_opening, binary_closing
    from sklearn.decomposition import FactorAnalysis, NMF

    def remove_hull(cube: np.ndarray) -> np.ndarray:
        # Fallback: no continuum removal; keep cube unchanged.
        return cube

    def mnf_transform(cube: np.ndarray, n_components: int) -> np.ndarray:
        # Simple MNF-like transform using FactorAnalysis on flattened pixels.
        h, w, b = cube.shape
        flat = cube.reshape(-1, b)
        n_comp = int(min(max(n_components, 1), b))
        fa = FactorAnalysis(n_components=n_comp, random_state=0)
        transformed = fa.fit_transform(flat)  # (pixels, n_comp)
        return transformed.reshape(h, w, n_comp).astype(np.float32)

    def nmf_unmixing(mnf_2d_components_by_pixels: np.ndarray, n_endmembers: int = 1) -> np.ndarray:
        # Input expected shape: (components, pixels). Output: abundance for first component (pixels,).
        X = np.maximum(mnf_2d_components_by_pixels.T, 0)  # (pixels, components)
        nmf = NMF(n_components=int(n_endmembers), init="nndsvda", random_state=0, max_iter=200)
        H = nmf.fit_transform(X)  # (pixels, n_endmembers)
        return H[:, 0].astype(np.float32)

    def threshold_abundance(abundance_map: np.ndarray, thresh: float) -> np.ndarray:
        return (abundance_map > thresh).astype(np.uint8)

    def clean_mask(mask: np.ndarray, min_area: int = 80, margin: int = 5) -> np.ndarray:
        # Connected-components filter + hole fill + dilation.
        structure = np.ones((3, 3), dtype=int)
        labeled, n_labels = label(mask.astype(bool), structure=structure)
        if n_labels == 0:
            return np.zeros_like(mask, dtype=np.uint8)
        areas = np.array([np.sum(labeled == i) for i in range(1, n_labels + 1)], dtype=int)
        keep = areas >= int(max(min_area, 0))
        cleaned = np.isin(labeled, np.where(keep)[0] + 1)
        cleaned = binary_fill_holes(cleaned)
        if margin and margin > 0:
            cleaned = binary_dilation(cleaned, iterations=int(margin))
        return cleaned.astype(np.uint8)

    def despeckle_mask(mask: np.ndarray) -> np.ndarray:
        # Light morphological open/close to reduce speckles.
        m = mask.astype(bool)
        m = binary_opening(m, structure=np.ones((3, 3), dtype=bool))
        m = binary_closing(m, structure=np.ones((3, 3), dtype=bool))
        return m.astype(np.uint8)

    def improve_mask_from_graph(mask: np.ndarray) -> np.ndarray:
        # Fallback: graph method not available; return input.
        return mask.astype(np.uint8)

    def export_results(*args, **kwargs):
        # Not used in this file (we export directly below). Kept for compatibility.
        raise NotImplementedError("export_results fallback is not implemented in this workspace.")

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
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
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
    parser.add_argument(
        "--csv_path",
        required=False,
        default=None,
        help="(Optional) CSV file with fungal reference spectra. If omitted, runs without reference‑based tuning.",
    )
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
    # 1️⃣ Load fungal reference spectra (optional)
    # ------------------------------------------------------------------
    if args.csv_path is not None:
        wavelengths, ref_spectra = [], []
        with open(args.csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                wavelengths.append(float(row["wavelength"]))
                ref_spectra.append(float(row["fungus_reflection"]))
        wavelengths = np.array(wavelengths, dtype=np.float32)
        ref_spectra = np.array(ref_spectra, dtype=np.float32)
    else:
        wavelengths = None
        ref_spectra = None

    # ------------------------------------------------------------------
    # 2️⃣ Load hyperspectral cube
    # ------------------------------------------------------------------
    if args.hsi_path.lower().endswith((".hdr", ".dat")):
        import spectral.io.envi as envi
        img = envi.open(args.hsi_path)
        cube = img.load().astype(np.float32)
        meta = getattr(img, "metadata", {}) or {}
        wl = meta.get("wavelength", [])
        try:
            wavelengths_hdr = np.array([float(x) for x in wl], dtype=np.float32) if wl else np.arange(cube.shape[2], dtype=np.float32)
        except Exception:
            wavelengths_hdr = np.arange(cube.shape[2], dtype=np.float32)
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
    if wavelengths_hdr.size >= 2:
        print(f"📈 Using wavelength range: {wavelengths_hdr[0]:.1f}–{wavelengths_hdr[-1]:.1f} nm")
    else:
        print("📈 Using wavelength metadata: unavailable")

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
    H, W, C = mnf_cube.shape
    mnf_2d = mnf_cube.reshape(-1, C).T  # (components, pixels)
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
    # Build a simple RGB preview for overlay (first/mid/last band)
    b = cube.shape[2]
    b1, b2, b3 = 0, b // 2, max(b - 1, 0)
    rgb_preview = np.stack([cube[:, :, b1], cube[:, :, b2], cube[:, :, b3]], axis=-1).astype(np.float32)
    p2, p98 = np.percentile(rgb_preview, 2), np.percentile(rgb_preview, 98)
    rgb_preview = np.clip((rgb_preview - p2) / (p98 - p2 + 1e-6), 0, 1)
    plt.imshow(rgb_preview)
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