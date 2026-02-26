# Quick Start Guide: Segmenting Cheese hyperspectral images (HSI) for Fungus Defects

This folder contains a **minimal synthetic test** and a **Python script** that implements an enhanced CoreSpecViewer‑based pipeline for segmenting fungal defect regions in cheese hyperspectral cubes.

## 1. Prerequisites

1. **Python 3.13+** (the environment already uses a virtual environment at  
   `.../.pi/venv`).

2. **Required Python packages** (most are pre‑installed; if missing install with pip):
   - `numpy`
   - `scipy`
   - `scikit-learn`
   - `spectral`
   - `rasterio`
   - `matplotlib`
   - `tqdm`
   - `opencv-python-headless`   *(required for some mask cleaning operations)*
   - `gfit` *(internal package – ensure it is on the PYTHONPATH)*  

   Example (if needed):
   ```bash
   pip install --break-system-packages opencv-python-headless
   ```

3. **Input data**:
   - **HSI cube**: ENVI file (`.hdr` + `.dat`) or a NumPy array saved with `.npy`.
   - **Reference CSV**: Must contain two columns named `wavelength` and `fungus_reflection`.  
     Example line: `2200,0.5231`

## 2. Output

Running the script produces, in the `--output_dir` you specify:

- `*_mask.tif` – binary mask (GeoTIFF) highlighting detected fungus pixels.
- `*_overlay.png` – visual overlay of the mask on an RGB preview of the cube.
- `summary.csv` – (optional) one‑row CSV with total pixels, fungus pixels, and percentage.

## 3. Running the Script

```bash
python3 segment_cheese_fungus_enhanced.py \
    --hsi_path /path/to/your_hsi.npy \
    --csv_path /path/to/reference.csv \
    --output_dir /path/to/results \
    [--thresh 0.12] \
    [--mnf_components 15] \
    [--min_mask_area 80] \
    [--margin 5] \
    [--clean_method none|despeckle|graph] \
    [--snr_thresh 20.0] \
    [--snr_min_run 20] \
    [--log_summary]
```

### Common flags

| Flag | Description |
|------|-------------|
| `--thresh` | Abundance threshold (default `0.12`). Lower → more pixels flagged. |
| `--mnf_components` | Number of MNF components (default `15`). |
| `--clean_method` | Mask cleaning strategy: `none` (no extra clean), `despeckle`, or `graph`. |
| `--log_summary` | Write a `summary.csv` with overall statistics. |
| `--snr_thresh` / `--snr_min_run` | Tune automatic band selection based on Signal‑to‑Noise Ratio. |

## 4. Tips for Successful Segmentation

1. **Band selection** – Fungal defects often have characteristic reflectance around **2200‑2320 nm**.  
   If you notice poor detection, try adjusting `--snr_thresh` or provide a custom wavelength range by editing the script.

2. **Threshold tuning** – Start with `--thresh 0.05` for low‑contrast defects, then increase if you get many false positives.

3. **MNF components** – 5‑10 components usually capture most variance; increase only if the cube is very noisy.

4. **Mask cleaning** –  
   - `despeckle` uses morphological opening to remove speckles.  
   - `graph` uses the `improve_mask_from_graph` heuristic for column‑wise thickening.

5. **Output inspection** – The generated overlay PNG is quick to glance at; the GeoTIFF mask can be loaded in any GIS viewer for precise area calculations.

## 5. Synthetic Test (already in this folder)

A tiny synthetic dataset was generated for verification:

- `synthetic.npy` – 20×20×8 cube with random reflectance.
- `synthetic_ref.csv` – reference spectra spanning 2200‑2300 nm.

Run the script on this data to confirm that the pipeline runs from start to finish:

```bash
python3 segment_cheese_fungus_enhanced.py \
    --hsi_path synthetic_test/synthetic.npy \
    --csv_path synthetic_test/synthetic_ref.csv \
    --output_dir synthetic_test/results \
    --thresh 0.05 \
    --mnf_components 5 \
    --clean_method none \
    --log_summary
```

If everything loads, you will see a mask, an overlay image, and a `summary.csv` printed in the console.

## 6. Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ModuleNotFoundError: No module named 'cv2'` | OpenCV not installed | `pip install --break-system-packages opencv-python-headless` |
| `ModuleNotFoundError: No module named 'gfit'` | Internal package missing from PYTHONPATH | Add the site‑packages directory to `PYTHONPATH` (see the script invocation example above). |
| Empty mask or all‑zero output | Threshold too high / insufficient MNF components | Lower `--thresh` or increase `--mnf_components`. |
| Import errors for `spectral` | Library not installed | `pip install spectral` |

---

**That's it!** You now have a ready‑to‑run script for detecting fungal defects in cheese hyperspectral images. Feel free to adapt the parameters to your specific dataset and enjoy the segmentation results. 🚀