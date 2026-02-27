#!/usr/bin/env python3
"""
Integration‑test script for cheese fungus segmentation.
* Reads a small synthetic HSI cube (test_cube.npy) and a reference CSV (test_fungus_ref.csv).
* Performs a very simple processing chain:
  - Minimum‑Noise‑Fraction via FactorAnalysis
  - NMF unmixing (single component)
  - Thresholding
  - Minimal cleaning using scipy.ndimage (binary dilation / fill_holes)
  - Export mask (GeoTIFF), overlay PNG, and summary CSV
* No external processing tools (e.g., CoreSpecViewer) are required,
  only numpy, pandas‑free csv, scipy, rasterio, matplotlib, tqdm.
"""

import os
import argparse
import csv
import numpy as np
import rasterio
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from scipy.ndimage import label, binary_fill_holes, binary_dilation, binary_erosion
from tqdm import tqdm
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import NMF

# ----------------------------------------------------------------------
#  Helper functions
# ----------------------------------------------------------------------
def load_reference(csv_path: str):
    """Read CSV with columns: wavelength, fungus_reflection."""
    wavelengths = []
    spectra = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            wavelengths.append(float(row["wavelength"]))
            spectra.append(float(row["fungus_reflection"]))
    return np.array(wavelengths, dtype=np.float32), np.array(spectra, dtype=np.float32)


def load_hsi(hsi_path: str):
    """
    Load hyperspectral cube.
    Supports:
    - ENVI files via .hdr/.dat (using spectral.io.envi)
    - .npy stacks (shape: rows, cols, bands)
    """
    lower = hsi_path.lower()
    if lower.endswith((".hdr", ".dat")):
        import spectral.io.envi as envi

        img = envi.open(hsi_path)
        data = img.load().astype(np.float32)
        # Wavelengths are not currently used downstream; provide a simple index array.
        wavelengths = np.arange(data.shape[2], dtype=np.float32)
    else:
        data = np.load(hsi_path).astype(np.float32)
        wavelengths = np.arange(data.shape[2], dtype=np.float32)
    return data, wavelengths


def threshold_abundance(abundance: np.ndarray, thresh: float) -> np.ndarray:
    """Binary mask where abundance exceeds thresh."""
    return (abundance > thresh).astype(np.uint8)


def clean_mask(mask: np.ndarray, min_area: int, margin: int):
    """
    Minimal cleaning:
    1. Label connected components (8‑connectivity).
    2. Keep components with size >= min_area.
    3. Fill holes.
    4. Optionally dilate by ``margin`` pixels.
    Returns a binary mask (0/1).
    """
    # 8‑connectivity structuring element (3×3 full ones)
    structure = np.ones((3, 3), dtype=int)

    # Label components
    labeled, n_labels = label(mask, structure=structure)

    # Compute area of each label (label 0 is background)
    areas = np.array([np.sum(labeled == i) for i in range(1, n_labels + 1)])
    keep = areas >= min_area

    # Build cleaned mask
    cleaned = np.isin(labeled, np.where(keep)[0] + 1)

    # Fill holes inside kept regions
    cleaned = binary_fill_holes(cleaned)

    # Dilate to add margin (if margin > 0)
    if margin > 0:
        cleaned = binary_dilation(cleaned, iterations=margin)

    return cleaned.astype(np.uint8)


def export_results(mask: np.ndarray, output_dir: str, base_name: str,
                  thresh: float, min_area: int, margin: int):
    """Write mask GeoTIFF, overlay PNG, and summary CSV."""
    mask_path = os.path.join(output_dir, f"{base_name}_mask.tif")
    overlay_path = os.path.join(output_dir, f"{base_name}_overlay.png")
    summary_path = os.path.join(output_dir, f"{base_name}_summary.csv")

    # ---- Mask GeoTIFF ----------------------------------------------------
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
    # Show a blank canvas; in a real run you could plot the RGB data here.
    plt.imshow(np.zeros((mask.shape[0], mask.shape[1], 3), dtype=np.uint8))
    # Draw mask in semi‑transparent red
    mask_vis = np.ma.masked_where(mask == 0, mask)
    plt.imshow(mask_vis, cmap="Reds", alpha=0.4)
    plt.tight_layout()
    plt.savefig(overlay_path, dpi=300, bbox_inches="tight")
    plt.close()

    # ---- Summary CSV ------------------------------------------------------
    total = mask.size
    fungus = mask.sum()
    percent = fungus / total * 100
    summary = {"total_pixels": total, "fungus_pixels": fungus,
               "percentage_fungus": percent}
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["total_pixels",
                                              "fungus_pixels",
                                              "percentage_fungus"])
        writer.writeheader()
        writer.writerow(summary)

    print(f"\nMask saved          -> {mask_path}")
    print(f"Overlay saved       -> {overlay_path}")
    print(f"Summary CSV saved   -> {summary_path}\n")


# ----------------------------------------------------------------------
#  Main routine
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Very simple fungus segmentation for testing.")
    parser.add_argument("--hsi_path", required=True,
                        help="Path to .npy hyperspectral cube.")
    parser.add_argument(
        "--csv_path",
        required=False,
        default=None,
        help="(Optional) Path to reference CSV (wavelength, fungus_reflection). Not required for this integration pipeline.",
    )
    parser.add_argument("--output_dir", default="./output",
                        help="Directory for output files.")
    parser.add_argument("--thresh", type=float, default=0.12,
                        help="Abundance threshold.")
    parser.add_argument("--mnf_components", type=int, default=5,
                        help="Number of FactorAnalysis components (MNF).")
    parser.add_argument("--min_mask_area", type=int, default=5,
                        help="Minimum component size (pixels).")
    parser.add_argument("--margin", type=int, default=2,
                        help="Margin (pixels) for dilation.")
    parser.add_argument("--roi_mask", default=None,
                        help="Optional .npy boolean mask (rows x cols). Only pixels with True are kept in final mask.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------- Load reference (optional) -----------------------
    if args.csv_path is not None:
        wavelengths, ref_spec = load_reference(args.csv_path)
        print(f"Reference spectra loaded: {len(wavelengths)} bands")
    else:
        print("Reference spectra: skipped (no --csv_path provided)")

    # ------------------- Load cube ---------------------------------------
    cube, band_centers = load_hsi(args.hsi_path)
    base_name = os.path.splitext(os.path.basename(args.hsi_path))[0]
    print(f"Cube loaded: {cube.shape[0]}x{cube.shape[1]}x{cube.shape[2]}")

    # ------------------- Minimum Noise Fraction ------------------------
    print(f"Computing {args.mnf_components} MNF components ...")
    # Reshape to (pixels, bands) for FactorAnalysis
    h, w, b = cube.shape
    flat = cube.reshape(-1, b)  # (pixels, bands)
    n_comp = min(args.mnf_components, b)
    fa = FactorAnalysis(n_components=n_comp, random_state=0)
    transformed = fa.fit_transform(flat)        # (pixels, n_comp)
    mnf = transformed.reshape(h, w, n_comp)     # (rows, cols, comps)

    # ------------------- NMF unmixing (single component) ---------------
    print("NMF unmixing ...")
    # Use the MNF components as features (pixels x components)
    pixels = h * w
    comps = mnf.shape[2]
    mnf_flat = mnf.reshape(pixels, comps)
    nmf = NMF(n_components=1, init="nndsvda", random_state=0, max_iter=200)
    H = nmf.fit_transform(np.maximum(mnf_flat, 0))  # (pixels, 1)
    abundance = H[:, 0].reshape(h, w)

    # ------------------- Threshold ---------------------------------------
    print(f"Thresholding (> {args.thresh}) ...")
    raw_mask = threshold_abundance(abundance, args.thresh)

    # ------------------- Clean mask --------------------------------------
    print("Cleaning mask ...")
    clean_mask_result = clean_mask(raw_mask, args.min_mask_area, args.margin)

    # ------------------- Apply ROI mask if provided ----------------------
    if args.roi_mask is not None:
        roi = np.load(args.roi_mask).astype(bool)
        if roi.shape != clean_mask_result.shape:
            raise ValueError(
                f"ROI mask shape {roi.shape} does not match cube spatial shape {clean_mask_result.shape}"
            )
        clean_mask_result = (clean_mask_result.astype(bool) & roi).astype(np.uint8)

    # ------------------- Export -------------------------------------------
    export_results(clean_mask_result, args.output_dir, base_name,
                   args.thresh, args.min_mask_area, args.margin)


if __name__ == "__main__":
    main()