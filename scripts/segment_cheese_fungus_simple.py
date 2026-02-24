#!/usr/bin/env python3
"""
Simple fungus segmentation script for cheese hyperspectral images.
Works with:
- hyperspectral cubes saved as .npy (numpy arrays) with shape (rows, cols, bands)
- reference CSV with columns: wavelength, fungus_reflection
Uses:
- Minimum Noise Fraction (FactorAnalysis)
- Non‑negative Matrix Factorization (NMF) to unmix the first end‑member
- Thresholding + morphological cleaning
- Exports mask (GeoTIFF), overlay PNG, and summary CSV
"""

import os
import argparse
import numpy as np
import csv
import rasterio
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from skimage import morphology, measure, filters
from skimage.morphology import disk
from tqdm import tqdm
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import NMF

# ----------------------------------------------------------------------
#  Helper functions
# ----------------------------------------------------------------------
def load_reference(csv_path: str):
    """Read reference CSV – returns wavelengths and spectra as np.ndarray."""
    wavelengths = []
    spectra = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            wavelengths.append(float(row["wavelength"]))
            spectra.append(float(row["fungus_reflection"]))
    return np.array(wavelengths, dtype=np.float32), np.array(spectra, dtype=np.float32)


def load_hsi(hsi_path: str):
    """Load a hyperspectral cube saved as .npy (rows, cols, bands)."""
    data = np.load(hsi_path).astype(np.float32)
    # No wavelength metadata – we just use the order of bands
    wavelengths = np.arange(data.shape[2], dtype=np.float32)
    return data, wavelengths


def threshold_abundance(abundance: np.ndarray, thresh: float) -> np.ndarray:
    """Binary mask where abundance > thresh."""
    return (abundance > thresh).astype(np.uint8)


def clean_mask(mask: np.ndarray, min_area: int, margin: int):
    """Remove small components, fill holes, optionally expand margin."""
    labeled = measure.label(mask)
    props = measure.regionprops(labeled)
    keep_labels = [p.label for p in props if p.area >= min_area]
    keep_mask = np.isin(labeled, keep_labels)

    filled = morphology.remove_small_holes(keep_mask, area_threshold=min_area, connectivity=2)

    if margin > 0:
        struct = disk(margin)
        filled = morphology.dilation(filled, struct)

    return filled.astype(np.uint8)


def export_results(mask: np.ndarray, output_dir: str, hsi_basename: str,
                  thresh: float, min_area: int, margin: int):
    """Write mask GeoTIFF, overlay PNG, and summary CSV."""
    mask_path = os.path.join(output_dir, f"{hsi_basename}_mask.tif")
    overlay_path = os.path.join(output_dir, f"{hsi_basename}_overlay.png")
    summary_path = os.path.join(output_dir, f"{hsi_basename}_summary.csv")

    # ---- Mask as GeoTIFF -------------------------------------------------
    meta = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 1,
        "height": mask.shape[0],
        "width": mask.shape[1],
        "transform": from_origin(0, 0, 1, 1)
    }
    os.makedirs(output_dir, exist_ok=True)
    with rasterio.open(mask_path, "w", **meta) as dst:
        dst.write(mask.astype(np.uint8), 1)

    # ---- Overlay PNG -----------------------------------------------------
    plt.figure(figsize=(8, 8))
    plt.axis("off")
    # just show a gray square; we don't have actual RGB data here
    plt.imshow(np.zeros((mask.shape[0], mask.shape[1], 3), dtype=np.uint8))
    # draw mask in semi‑transparent red
    mask_vis = np.ma.masked_where(mask == 0, mask)
    plt.imshow(mask_vis, cmap="Reds", alpha=0.4)
    plt.tight_layout()
    plt.savefig(overlay_path, dpi=300, bbox_inches="tight")
    plt.close()

    # ---- Summary CSV ------------------------------------------------------
    total_pixels = mask.size
    fungus_pixels = mask.sum()
    percent = fungus_pixels / total_pixels * 100
    summary = {
        "total_pixels": total_pixels,
        "fungus_pixels": fungus_pixels,
        "percentage_fungus": percent,
    }
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["total_pixels", "fungus_pixels", "percentage_fungus"])
        writer.writeheader()
        writer.writerow(summary)

    print(f"\n✅ Mask saved          → {mask_path}")
    print(f"✅ Overlay saved       → {overlay_path}")
    print(f"✅ Summary CSV saved   → {summary_path}\n")


# ----------------------------------------------------------------------
#  Main routine
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Simple fungus defect segmentation in cheese HSI."
    )
    parser.add_argument("--hsi_path", required=True,
                        help="Path to hyperspectral cube saved as .npy.")
    parser.add_argument("--csv_path", required=True,
                        help="CSV file with fungal reference spectra.")
    parser.add_argument("--output_dir", default="./results",
                        help="Folder for output files.")
    parser.add_argument("--thresh", type=float, default=0.12,
                        help="Abundance threshold for mask.")
    parser.add_argument("--mnf_components", type=int, default=5,
                        help="Number of MNF components (FactorAnalysis).")
    parser.add_argument("--min_mask_area", type=int, default=5,
                        help="Minimum object size (pixels) after cleaning.")
    parser.add_argument("--margin", type=int, default=2,
                        help="Margin (pixels) for overlay expansion.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------- Load reference spectra -------------------------
    wavelengths, ref_spectra = load_reference(args.csv_path)
    print(f"🔍 Loaded reference: {len(wavelengths)} bands")

    # ------------------- Load hyperspectral cube -----------------------
    cube, band_centers = load_hsi(args.hsi_path)
    hsi_basename = os.path.splitext(os.path.basename(args.hsi_path))[0]
    print(f"🔍 Loaded cube: {cube.shape[0]}×{cube.shape[1]}×{cube.shape[2]}")

    # ------------------- Minimum Noise Fraction -----------------------
    print(f"⚙️  Computing {args.mnf_components} MNF components …")
    # Center the data (FactorAnalysis assumes zero‑mean)
    centered = cube.reshape(-1, cube.shape[2]).T  # (bands, pixels)
    fa = FactorAnalysis(n_components=args.mnf_components, random_state=0)
    fa.fit(centered)
    transformed = fa.transform(centered)  # (components, pixels)
    mnf = transformed.T.reshape(cube.shape)  # (rows, cols, comps)
    mnf_2d = mnf.reshape(-1, args.mnf_components).T  # (components, pixels)

    # ------------------- NMF unmixing (first component) ---------------
    print("🧩  NMF unmixing …")
    nmf = NMF(n_components=1, random_state=0, alpha=0.1)
    # Fit on the transpose because NMF expects (features, samples)
    W = nmf.fit_transform(mnf_2d.T)   # (1, pixels)
    H = nmf.transform(mnf_2d)         # (pixels, 1)
    abundance_map = H[:, 0].reshape(cube.shape[:2])

    # ------------------- Threshold → binary mask ----------------------
    print(f"🔎  Thresholding (>{args.thresh}) …")
    raw_mask = threshold_abundance(abundance_map, args.thresh)

    # ------------------- Clean mask ------------------------------------
    print("🧹  Cleaning mask …")
    clean_mask_result = clean_mask(raw_mask, min_area=args.min_mask_area, margin=args.margin)

    # ------------------- Export results --------------------------------
    export_results(clean_mask_result, args.output_dir, hsi_basename,
                   args.thresh, args.min_mask_area, args.margin)

    print("✅ All done!")


if __name__ == "__main__":
    main()