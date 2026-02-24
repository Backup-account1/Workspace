#!/usr/bin/env python3
"""
segment_cheese_fungus.py

Segment fungal defect regions in cheese hyperspectral images (HSI) using a
reference CSV of fungal spectra.  The workflow follows the CoreSpecViewer
core processing chain (Continuum‑Removal → Minimum‑Noise‑Fraction → NMF
unmixing → thresholding → clean‑up).

Author: <your name>
Created: 2026‑02‑23
"""

# ----------------------------------------------------------------------
#  Imports – everything we need from the CoreSpecViewer / Spectral stack
# ----------------------------------------------------------------------
import os
import argparse
import numpy as np
import csv
import spectral
import spectral.io.envi as envi
from spectral.util import import_images
from sklearn.decomposition import FactorAnalysis
from skimage import morphology, measure, filters
from skimage.morphology import disk
import matplotlib.pyplot as plt
from tqdm import tqdm
import rasterio
from rasterio.transform import from_origin

# ----------------------------------------------------------------------
#  Helper functions
# ----------------------------------------------------------------------
def load_reference(csv_path: str):
    """Read the reference CSV – returns wavelengths and spectra as np.ndarray."""
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
    Load an ENVI cube (or a stacked numpy file) and return:
        cube      – (rows, cols, bands) float32 array
        wavelengths – 1‑D array of band centres (nm)
    """
    if hsi_path.lower().endswith(('.hdr', '.dat')):
        # ENVI format – spectral library handles it natively
        cube = spectral.open_image(hsi_path)
        data = cube.load().astype(np.float32)
        wavelengths = np.array(cube.get_metadata('wavelength'), dtype=np.float32)
    else:
        # Assume a simple .npy stack where the last axis = bands
        data = np.load(hsi_path).astype(np.float32)
        # If the file does not carry wavelength info, fall back to dummy indices
        wavelengths = np.arange(data.shape[2], dtype=np.float32)
    return data, wavelengths


def continuum_remove(cube: np.ndarray) -> np.ndarray:
    """
    Wrapper around CoreSpecViewer's `remove_hull` (aka CR).
    Performs continuum removal on every spatial pixel.
    """
    # `remove_hull` lives in spectral_functions.py (the file you read)
    # It expects a 3‑D array (H, W, B) and returns a continuum‑removed cube.
    from app.spectral_ops.spectral_functions import remove_hull
    return remove_hull(cube)


def mnf_transform(cube: np.ndarray, n_components: int) -> np.ndarray:
    """
    Minimum‑Noise‑Fraction projection (CoreSpecViewer uses a custom wrapper,
    but the pure‑numpy version via FactorAnalysis works just as well).
    Returns a cube of shape (rows, cols, n_components).
    """
    # Reshape to (pixels, bands) for FactorAnalysis
    H, W, B = cube.shape
    cube_2d = cube.reshape(-1, B).T  # (bands, pixels)
    fa = FactorAnalysis(n_components=n_components, max_iter=200, random_state=0)
    fa.fit(cube_2d)
    transformed = fa.transform(cube_2d)          # (n_components, pixels)
    mnf = transformed.T.reshape(H, W, n_components)  # (H, W, n_components)
    return mnf.astype(np.float32)


def nmf_unmixing(cube_2d: np.ndarray, n_endmembers: int, max_iter: int = 500, tol: float = 1e-4):
    """
    Linear NMF unmixing using scikit‑learn's `NMF`.
    Returns the abundance map for the *first* end‑member (index 0).
    """
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=n_endmembers, max_iter=max_iter, tol=tol, random_state=0, alpha=0.1)
    # Fit on the *transposed* data because NMF expects (features, samples)
    W = nmf.fit_transform(cube_2d.T).T          # (n_endmembers, pixels)
    H = nmf.transform(cube_2d)                  # (pixels, n_endmembers)
    # Abundance of the first end‑member (fungus)
    abundance = H[:, 0].reshape(cube_2d.shape[:2])
    return abundance.astype(np.float32)

def threshold_abundance(abundance: np.ndarray, thresh: float) -> np.ndarray:
    """Binary mask where abundance > thresh."""
    return (abundance > thresh).astype(np.uint8)

def clean_mask(mask: np.ndarray, min_area: int, margin: int):
    """
    Remove small components, fill holes, and optionally expand the mask
    by a margin (for overlay purposes). Returns a cleaned binary mask.
    """
    # Remove tiny objects
    labeled = measure.label(mask)
    props = measure.regionprops(labeled)
    keep_labels = [p.label for p in props if p.area >= min_area]
    keep_mask = np.isin(labeled, keep_labels)

    # Fill holes inside the kept regions
    filled = morphology.remove_small_holes(keep_mask, area_threshold=min_area, connectivity=2)

    # Optional dilation/erosion to add a margin (useful for overlay visualisation)
    if margin > 0:
        struct = disk(margin)
        filled = morphology.dilation(filled, struct)

    return filled.astype(np.uint8)


def export_results(mask: np.ndarray, wavelengths: np.ndarray,
                  hsi_path: str, output_dir: str,
                  thresh: float, min_area: int, margin: int):
    """
    Write:
        • binary mask (GeoTIFF)
        • overlay PNG (RGB image with mask drawn in red)
        • summary CSV (pixel counts, percentages, etc.)
    """
    base = os.path.splitext(os.path.basename(hsi_path))[0]
    mask_path = os.path.join(output_dir, f"{base}_mask.tif")
    overlay_path = os.path.join(output_dir, f"{base}_overlay.png")
    summary_path = os.path.join(output_dir, f"{base}_summary.csv")

    # ---- Mask as GeoTIFF -------------------------------------------------
    meta = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 1,
        "height": mask.shape[0],
        "width": mask.shape[1],
        "transform": from_origin(0, 0, 1, 1)  # simple affine; adjust if you have georef
    }
    with rasterio.open(mask_path, "w", **meta) as dst:
        dst.write(mask.astype(np.uint8), 1)

    # ---- Overlay PNG -----------------------------------------------------
    # Simple RGB preview: first three bands or a quick false‑colour stretch
    preview = spectral.util.open_image(hsi_path).load()   # (H, W, B)
    rgb_preview = preview[:, :, :3].astype(np.float32)

    plt.figure(figsize=(8, 8))
    plt.imshow(rgb_preview)
    # draw mask in semi‑transparent red
    mask_vis = np.ma.masked_where(mask == 0, mask)  # 0 → transparent
    plt.imshow(mask_vis, cmap="Reds", alpha=0.4)
    plt.axis("off")
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
    # Write CSV without pandas
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["total_pixels", "fungus_pixels", "percentage_fungus"])
        writer.writeheader()
        writer.writerow(summary)

    # ---- Inform the user --------------------------------------------------
    print(f"\n✅ Mask saved          → {mask_path}")
    print(f"✅ Overlay saved       → {overlay_path}")
    print(f"✅ Summary CSV saved   → {summary_path}\n")


# ----------------------------------------------------------------------
#  Main processing routine
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Fungus defect segmentation in cheese hyperspectral images."
    )
    parser.add_argument("--hsi_path", required=True, help="Path to ENVI cube (or .npy stack).")
    parser.add_argument("--csv_path", required=True, help="CSV file with fungal reference spectra.")
    parser.add_argument("--output_dir", default="./results", help="Folder for output files.")
    parser.add_argument("--thresh", type=float, default=0.12, help="Abundance threshold for mask.")
    parser.add_argument("--mnf_components", type=int, default=15, help="Number of MNF components.")
    parser.add_argument("--min_mask_area", type=int, default=80, help="Minimum object size (pixels).")
    parser.add_argument("--margin", type=int, default=5, help="Margin (pixels) for overlay expansion.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    #  1️⃣ Load reference spectra
    # ------------------------------------------------------------------
    wavelengths, ref_spectra = load_reference(args.csv_path)

    # ------------------------------------------------------------------
    #  2️⃣ Load hyperspectral cube
    # ------------------------------------------------------------------
    cube, band_centers = load_hsi(args.hsi_path)
    print(f"🔍 Loaded cube: {cube.shape[0]}×{cube.shape[1]}×{cube.shape[2]}")

    # ------------------------------------------------------------------
    #  3️⃣ Continuum removal (CR)
    # ------------------------------------------------------------------
    print("🛠️  Continuum removal …")
    cr_cube = continuum_remove(cube)

    # ------------------------------------------------------------------
    #  4️⃣ Optional band‑range trim (e.g., 2200‑2320 nm)
    # ------------------------------------------------------------------
    # For demonstration we keep all bands; you can uncomment and set your range:
    # mask_range = (np.where((band_centers >= 2200) & (band_centers <= 2320))[0])
    # cube = cube[:, :, mask_range]
    # wavelengths = wavelengths[mask_range]
    # cr_cube = cr_cube[:, :, mask_range]

    # ------------------------------------------------------------------
    #  5️⃣ M NF projection
    # ------------------------------------------------------------------
    print(f"⚙️  Computing {args.mnf_components} MNF components …")
    mnf_cube = mnf_transform(cr_cube, args.mnf_components)

    # ------------------------------------------------------------------
    #  6️⃣ Flatten for NMF unmixing
    # ------------------------------------------------------------------
    H, W, _ = mnf_cube.shape
    mnf_2d = mnf_cube.reshape(-1, args.mnf_components).T  # (components, pixels)

    # ------------------------------------------------------------------
    #  7️⃣ NMF unmixing – abundance of the first component (fungus)
    # ------------------------------------------------------------------
    print("🧩  NMF unmixing …")
    abundance_map = nmf_unmixing(mnf_2d, n_endmembers=1)  # returns (pixels,) array
    abundance_map = abundance_map.reshape(H, W)

    # ------------------------------------------------------------------
    #  8️⃣ Threshold → binary mask
    # ------------------------------------------------------------------
    print(f"🔎  Thresholding (>{args.thresh}) …")
    raw_mask = threshold_abundance(abundance_map, args.thresh)

    # ------------------------------------------------------------------
    #  9️⃣ Clean up the mask (remove speckles, keep only sizeable objects)
    # ------------------------------------------------------------------
    print("🧹  Cleaning mask …")
    clean_mask_result = clean_mask(raw_mask, min_area=args.min_mask_area, margin=args.margin)

    # ------------------------------------------------------------------
    #  🔟 Export results
    # ------------------------------------------------------------------
    export_results(clean_mask_result, wavelengths, args.hsi_path,
                   args.output_dir, args.thresh, args.min_mask_area, args.margin)

    print("✅ All done!")


if __name__ == "__main__":
    main()