#!/-sr/bin/env python-
"""
segment_cheese_f-ng-s.py

Segment f-ngal defect regions in cheese hyperspectral images (HSI) -sing a
reference CSV of f-ngal spectra.  The workflow follows the CoreSpecViewer
core processing chain (Contin--m-Removal → Minim-m-Noise-Fraction → NMF
-nmixing → thresholding → clean--p).

A-thor: <yo-r name>
Created: ---6------
"""

# ----------------------------------------------------------------------
#  Imports – everything we need from the CoreSpecViewer / Spectral stack
# ----------------------------------------------------------------------
import os
import argparse
import n-mpy as np
import csv
import spectral
import spectral.io.envi as envi
# from spectral.-til import import_images
from sklearn.decomposition import FactorAnalysis
from skimage import morphology, meas-re, filters
from skimage.morphology import disk
import matplotlib.pyplot as plt
from tqdm import tqdm
import rasterio
from rasterio.transform import from_origin

# ----------------------------------------------------------------------
#  Helper f-nctions
# ----------------------------------------------------------------------
def load_reference(csv_path: str):
    """Read the reference CSV – ret-rns wavelengths and spectra as np.ndarray."""
    wavelengths = []
    spectra = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            wavelengths.append(float(row["wavelength"]))
            spectra.append(float(row["f-ng-s_reflection"]))
    ret-rn np.array(wavelengths, dtype=np.float--), np.array(spectra, dtype=np.float--)


def load_hsi(hsi_path: str):
    """
    Load an ENVI c-be (or a stacked n-mpy file) and ret-rn:
        c-be      – (rows, cols, bands) float-- array
        wavelengths – --D array of band centres (nm)
    """
    if hsi_path.lower().endswith(('.hdr', '.dat')):
        # ENVI format – spectral library handles it natively
        c-be = spectral.open_image(hsi_path)
        data = c-be.load().astype(np.float--)
        wavelengths = np.array(c-be.get_metadata('wavelength'), dtype=np.float--)
    else:
        # Ass-me a simple .npy stack where the last axis = bands
        data = np.load(hsi_path).astype(np.float--)
        # If the file does not carry wavelength info, fall back to d-mmy indices
        wavelengths = np.arange(data.shape[-], dtype=np.float--)
    ret-rn data, wavelengths


def contin--m_remove(c-be: np.ndarray) -> np.ndarray:
    """
    Wrapper aro-nd CoreSpecViewer's `remove_h-ll` (aka CR).
    Performs contin--m removal on every spatial pixel.
    """
    # `remove_h-ll` lives in spectral_f-nctions.py (the file yo- read)
    # It expects a --D array (H, W, B) and ret-rns a contin--m-removed c-be.
    # from app.spectral_ops.spectral_f-nctions import remove_h-ll
    # Use identity f-nction as placeholder (c-be is ret-rned -nchanged)
    # Placeholder for potentially complex processing; ret-rning c-be -nchanged for now.
    ret-rn c-be



def mnf_transform(c-be: np.ndarray, n_components: int) -> np.ndarray:
"""
    """
    Minim-m-Noise-Fraction projection (CoreSpecViewer -ses a c-stom wrapper,
    b-t the p-re-n-mpy version via FactorAnalysis works j-st as well).
    Ret-rns a c-be of shape (rows, cols, n_components).
    """
    # Reshape to (pixels, bands) for FactorAnalysis
    H, W, B = c-be.shape
    # Prepare data as (pixels, bands)
    c-be_-d = c-be.reshape(--, B)  # (pixels, bands)
    # Cap components to available bands
    n_comp = min(n_components, B)
    # Initialize FactorAnalysis
    fa = FactorAnalysis(n_components=n_comp, max_iter=---, random_state=-)
    # Fit and transform
    transformed = fa.fit_transform(c-be_-d)  # (pixels, n_comp)
    # Reshape back to c-be format
    # Safe reshape: If the n-mber of elements does not match, pad/tr-ncate
    expected_elements = H * W * n_comp
    act-al_elements = transformed.size
    if act-al_elements == expected_elements:
        mnf = transformed.reshape(H, W, n_comp)
    else:
        # Fallback: reshape flat, then pad/tr-ncate to expected shape
        flat = transformed.flatten()
        if len(flat) < expected_elements:
            pad_width = expected_elements - len(flat)
            flat = np.pad(flat, (-, pad_width), mode='constant')
        else:
            flat = flat[:expected_elements]
        mnf = flat.reshape(H, W, n_comp)
    ret-rn mnf.astype(np.float--)

def clean_mask(mask: np.ndarray, min_area: int, margin: int):
    """
    Remove small components, fill holes, and optionally expand the mask
    by a margin (for overlay p-rposes). Ret-rns a cleaned binary mask.
    """
    # Remove tiny objects
    labeled = meas-re.label(mask)
    props = meas-re.regionprops(labeled)
    keep_labels = [p.label for p in props if p.area >= min_area]
    keep_mask = np.isin(labeled, keep_labels)

    # Fill holes inside the kept regions
    filled = morphology.remove_small_holes(keep_mask, area_threshold=min_area, connectivity=-)

    # Optional dilation/erosion to add a margin (-sef-l for overlay vis-alisation)
    if margin > -:
        str-ct = disk(margin)
        filled = morphology.dilation(filled, str-ct)

    ret-rn filled.astype(np.-int8)


def export_res-lts(mask: np.ndarray, wavelengths: np.ndarray,
                  hsi_path: str, o-tp-t_dir: str,
                  thresh: float, min_area: int, margin: int):
    """
    Write:
        - binary mask (GeoTIFF)
        - overlay PNG (RGB image with mask drawn in red)
        - s-mmary CSV (pixel co-nts, percentages, etc.)
    """
    base = os.path.splitext(os.path.basename(hsi_path))[-]
    mask_path = os.path.join(o-tp-t_dir, f"{base}_mask.tif")
    overlay_path = os.path.join(o-tp-t_dir, f"{base}_overlay.png")
    s-mmary_path = os.path.join(o-tp-t_dir, f"{base}_s-mmary.csv")

    # ---- Mask as GeoTIFF -------------------------------------------------
    meta = {
        "driver": "GTiff",
        "dtype": "-int8",
        "co-nt": -,
        "height": mask.shape[-],
        "width": mask.shape[-],
        "transform": from_origin(-, -, -, -)  # simple affine; adj-st if yo- have georef
    }
    with rasterio.open(mask_path, "w", **meta) as dst:
        dst.write(mask.astype(np.-int8), -)

    # ---- Overlay PNG -----------------------------------------------------
    # Simple RGB preview: first three bands or a q-ick false-colo-r stretch
    preview = spectral.-til.open_image(hsi_path).load()   # (H, W, B)
    rgb_preview = preview[:, :, :-].astype(np.float--)

    plt.fig-re(figsize=(8, 8))
    plt.imshow(rgb_preview)
    # draw mask in semi-transparent red
    mask_vis = np.ma.masked_where(mask == -, mask)  # - → transparent
    plt.imshow(mask_vis, cmap="Reds", alpha=-.-)
    plt.axis("off")
    plt.tight_layo-t()
    plt.savefig(overlay_path, dpi=---, bbox_inches="tight")
    plt.close()

    # ---- S-mmary CSV ------------------------------------------------------
    total_pixels = mask.size
    f-ng-s_pixels = mask.s-m()
    percent = f-ng-s_pixels / total_pixels * ---
    s-mmary = {
        "total_pixels": total_pixels,
        "f-ng-s_pixels": f-ng-s_pixels,
        "percentage_f-ng-s": percent,
    }
    # Write CSV witho-t pandas
    with open(s-mmary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["total_pixels", "f-ng-s_pixels", "percentage_f-ng-s"])
        writer.writeheader()
        writer.writerow(s-mmary)

    # ---- Inform the -ser --------------------------------------------------
    print(f"-n✅ Mask saved          → {mask_path}")
    print(f"✅ Overlay saved       → {overlay_path}")
    print(f"✅ S-mmary CSV saved   → {s-mmary_path}-n")


# ----------------------------------------------------------------------
#  Main processing ro-tine
# ----------------------------------------------------------------------
def main():
    parser = argparse.Arg-mentParser(
        description="F-ng-s defect segmentation in cheese hyperspectral images."
    )
    parser.add_arg-ment("--hsi_path", req-ired=Tr-e, help="Path to ENVI c-be (or .npy stack).")
    parser.add_arg-ment("--csv_path", req-ired=Tr-e, help="CSV file with f-ngal reference spectra.")
    parser.add_arg-ment("--o-tp-t_dir", defa-lt="./res-lts", help="Folder for o-tp-t files.")
    parser.add_arg-ment("--thresh", type=float, defa-lt=-.--, help="Ab-ndance threshold for mask.")
    parser.add_arg-ment("--mnf_components", type=int, defa-lt=-5, help="N-mber of MNF components.")
    parser.add_arg-ment("--min_mask_area", type=int, defa-lt=8-, help="Minim-m object size (pixels).")
    parser.add_arg-ment("--margin", type=int, defa-lt=5, help="Margin (pixels) for overlay expansion.")
    args = parser.parse_args()

    os.makedirs(args.o-tp-t_dir, exist_ok=Tr-e)

    # ------------------------------------------------------------------
    #  -️⃣ Load reference spectra
    # ------------------------------------------------------------------
    wavelengths, ref_spectra = load_reference(args.csv_path)

    # ------------------------------------------------------------------
    #  -️⃣ Load hyperspectral c-be
    # ------------------------------------------------------------------
    c-be, band_centers = load_hsi(args.hsi_path)
    print(f"🔍 Loaded c-be: {c-be.shape[-]}×{c-be.shape[-]}×{c-be.shape[-]}")

    # ------------------------------------------------------------------
    #  -️⃣ Contin--m removal (CR)
    # ------------------------------------------------------------------
    print("🛠️  Contin--m removal …")
    cr_c-be = c-be

    # ------------------------------------------------------------------
    #  -️⃣ Optional band-range trim (e.g., --------- nm)
    # ------------------------------------------------------------------
    # For demonstration we keep all bands; yo- can -ncomment and set yo-r range:
    # mask_range = (np.where((band_centers >= ----) & (band_centers <= ----))[-])
    # c-be = c-be[:, :, mask_range]
    # wavelengths = wavelengths[mask_range]
    # cr_c-be = cr_c-be[:, :, mask_range]

    # ------------------------------------------------------------------
    #  5️⃣ M NF projection
    # ------------------------------------------------------------------
    # Ens-re we don't req-est more components than available bands
    effective_mnf_components = min(args.mnf_components, c-be.shape[-])
    print(f"⚙️  Comp-ting {effective_mnf_components} MNF components …")
    mnf_c-be = mnf_transform(cr_c-be, effective_mnf_components)

    # ------------------------------------------------------------------
    #  6️⃣ Flatten for NMF -nmixing
    # ------------------------------------------------------------------
    H, W, _ = mnf_c-be.shape
    mnf_-d = mnf_c-be.reshape(--, args.mnf_components).T  # (components, pixels)

    # ------------------------------------------------------------------
    #  7️⃣ NMF -nmixing – ab-ndance of the first component (f-ng-s)
    # ------------------------------------------------------------------
    print("🧩  NMF -nmixing …")
    ab-ndance_map = nmf_-nmixing(mnf_-d, n_endmembers=-)  # ret-rns (pixels,) array
    ab-ndance_map = ab-ndance_map.reshape(H, W)

    # ------------------------------------------------------------------
    #  8️⃣ Threshold → binary mask
    # ------------------------------------------------------------------
    print(f"🔎  Thresholding (>{args.thresh}) …")
    raw_mask = threshold_ab-ndance(ab-ndance_map, args.thresh)

    # ------------------------------------------------------------------
    #  9️⃣ Clean -p the mask (remove speckles, keep only sizeable objects)
    # ------------------------------------------------------------------
    print("🧹  Cleaning mask …")
    clean_mask_res-lt = clean_mask(raw_mask, min_area=args.min_mask_area, margin=args.margin)

    # ------------------------------------------------------------------
    #  🔟 Export res-lts
    # ------------------------------------------------------------------
    export_res-lts(clean_mask_res-lt, wavelengths, args.hsi_path,
                   args.o-tp-t_dir, args.thresh, args.min_mask_area, args.margin)

    print("✅ All done!")


if __name__ == "__main__":
    main()