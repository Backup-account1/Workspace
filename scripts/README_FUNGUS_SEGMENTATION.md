# 📚 README – Fungus Segmentation in Cheese Hyperspectral Images

**Tool:** `segment_cheese_fungus.py`   
**Location:** `~/openclaw/workspace/scripts/`   

This script automates the detection of fungal defects in cheese hyperspectral images (HSI) using a reference spectra CSV and the CoreSpecViewer processing chain (Continuum‑Removal → Minimum‑Noise‑Fraction → NMF unmixing → thresholding → morphological cleaning).  

---  

## Table of Contents
1. [Prerequisites](#prerequisites)  
2. [CSV Reference File Format](#csv-reference-file-format)  
3. [Installation & Environment Setup](#installation--environment-setup)  
4. [Running the Script](#running-the-script)  
5. [Parameter Details & Tips](#parameter-details--tips)  
6. [Output Files](#output-files)  
7. [Troubleshooting & FAQs](#troubleshooting--faqs)  
8. [Citation & Acknowledgments](#citation--acknowledgments)  

---  

## Prerequisites

| Requirement | Command |
|------------|---------|
| **Python 3.9+** (recommended 3.10‑3.12) | `python --version` |
| **OpenClaw virtual environment** (contains CoreSpecViewer & Spectral) | `source ~/.openclaw/workspace/.pi/venv/bin/activate` |
| **Core libraries** (NumPy, Pandas, scikit‑image, rasterio, tqdm, matplotlib) | `pip install numpy pandas scikit-image rasterio tqdm matplotlib` |
| **Spectral Python** (already in the OpenClaw env) | `pip show spectral` |
| **Optional:** `rasterio` for GeoTIFF export | included above |

> **Note:** All the above packages are pre‑installed in the OpenClaw environment, so you typically only need to activate it.

---  

## CSV Reference File Format

The reference CSV defines the fungal spectral signature you want to detect.

| Column | Description | Required? |
|--------|-------------|----------|
| `wavelength` | Band centre wavelength in **nanometres** (must be sorted ascending). | ✅ |
| `fungus_reflection` | Reflectance value for that band (0 – 1). One row per training pixel of the fungus. | ✅ |
| `label` *(optional)* | Human‑readable label – ignored by the script. | ❌ |

**Example (`fungus_ref.csv`)**

```csv
wavelength, fungus_reflection
400,0.018
410,0.021
420,0.024
...
800,0.067
810,0.070
820,0.074
```

- The number of rows **must** equal the total number of spectral bands in the HSI file.  
- If your HSI uses a different order, reorder the CSV accordingly or edit the script to match.

---  

## Installation & Environment Setup

```bash
# 1️⃣ Activate the OpenClaw virtual environment
source ~/.openclaw/workspace/.pi/venv/bin/activate

# 2️⃣ (Optional) Verify required packages are present
pip list | grep -E "numpy|pandas|scikit-image|rasterio|tqdm|matplotlib|spectral"

# 3️⃣ Navigate to the scripts folder (where the script and README live)
cd ~/openclaw/workspace/scripts
```

---  

## Running the Script

```bash
python segment_cheese_fungus.py \
    --hsi_path      path/to/cheese_cube.hdr \
    --csv_path      path/to/fungus_ref.csv \
    --output_dir    ./fungus_results \
    --thresh        0.12 \
    --mnf_components 15 \
    --min_mask_area 80 \
    --margin        5
```

### Mandatory arguments
| Flag | Meaning |
|------|----------|
| `--hsi_path` | Path to the hyperspectral cube (ENVI `.hdr`/`.dat` or a stacked `.npy`). |
| `--csv_path` | Path to the reference CSV file (see format above). |

### Optional arguments (defaults in brackets)

| Flag | Default | Description |
|------|---------|-------------|
| `--output_dir` | `./results` | Folder where mask, overlay PNG, and summary CSV will be written. |
| `--thresh` | `0.12` | Abundance‑threshold; increase to get stricter detection, decrease for more inclusive masks. |
| `--mnf_components` | `15` | Number of Minimum‑Noise‑Fraction components kept. More components keep more variance but increase runtime. |
| `--min_mask_area` | `80` | Minimum object size (in pixels) after connected‑component filtering; helps discard speckles. |
| `--margin` | `5` | Extra pixel margin added when creating the overlay PNG (useful for visual context). |

---  

## Parameter Details & Tips

| Parameter | How it influences results | Practical tips |
|-----------|--------------------------|----------------|
| `--thresh` | Controls the **abundance cut‑off** before binarization. Lower values → larger mask (more pixels labelled fungus) but also more false positives. Higher values → tighter detection, risking missed defects. | Start with `0.1`–`0.15`. Visualise the output and adjust until the mask covers the expected fungus blobs without noisy speckles. |
| `--mnf_components` | Determines how many **noise‑reduced spectral dimensions** are fed to NMF. Too few → insufficient separation of fungus from background; too many → noise may dominate. | 10‑20 is typical for cheese HSI (≈ 200‑400 bands). If you increase the number of bands dramatically, you may need more components. |
| `--min_mask_area` | Filters out **tiny disconnected components** after thresholding. Larger values clean up noisy speckles but may also erase small genuine defects. | Set to roughly the expected minimum defect size (in pixels). Use the overlay PNG to gauge size. |
| `--margin` | Adds extra border around the binary mask when generating the overlay image. | Helpful when you want the mask to be clearly visible on the RGB preview; does not affect the saved mask file. |

---  

## Output Files

All outputs are written to the folder specified by `--output_dir`.

| File | Description |
|------|-------------|
| `<basename>_mask.tif` | Binary mask (uint8) where **1** = fungus‑detected pixel, **0** = background. Can be loaded in GIS software or visualised with `rasterio`. |
| `<basename>_overlay.png` | RGB preview of the HSI with the mask drawn in semi‑transparent red. Ideal for quick visual inspection. |
| `<basename>_summary.csv` | CSV summary containing total pixel count, fungus pixel count, % fungus, and mean/median abundance values for statistics. |

`<basename>` is derived from the input HSI filename (e.g., `cheese_cube_mask.tif` for `cheese_cube.hdr`).

---  

## Troubleshooting & FAQs

| Issue | Likely cause | Fix |
|-------|--------------|-----|
| **Script aborts on `load_hsi`** | Wrong file extension or missing ENVI header. | Ensure the cube is either an ENVI header+data pair or a `.npy` file created with `np.save`. Verify file accessibility (`ls -l <path>`). |
| **Empty or all‑zero mask** | Threshold too high, or reference spectra do not match the fungus in the image. | Lower `--thresh` (e.g., `0.05`) and double‑check that the CSV wavelengths align with the HSI bands. |
| **Very noisy mask with many speckles** | MNF components insufficient or `min_mask_area` too low. | Increase `--mnf_components` or raise `--min_mask_area`. |
| **Overlay PNG looks blank** | The HSI preview extraction failed (e.g., insufficient bands). | Ensure the cube has at least 3 bands; otherwise modify the script to use a different band selection (`preview = preview[:, :, :3]` can be changed). |
| **Runtime too long (>5 min)** | Large data cube (e.g., >500 MB) or high `--mnf_components`. | Consider down‑sampling spatially (e.g., `cube = cube[::2, ::2, :]`) for quick tests, or increase `--mnf_components` only if needed. |
| **`remove_hull` ImportError** | Running outside the OpenClaw environment. | Make sure you activated the virtual environment (`source ~/.openclaw/workspace/.pi/venv/bin/activate`). |

---  

## Citation & Acknowledgments

- **CoreSpecViewer** – Russell Rogers, Geological Survey Ireland. Open‑source hyperspectral processing toolkit.  
- **Spectral Python (SPy)** – https://github.com/spectralpython/spectral  
- **OpenClaw** – Open‑source AI assistant framework (https://openclaw.ai)  

If you use this script for published work, please cite the CoreSpecViewer paper and acknowledge the OpenClaw project.

---  

### Quick Command Cheat‑Sheet

```bash
# Activate environment
source ~/.openclaw/workspace/.pi/venv/bin/activate

# Run detection (replace placeholders)
python segment_cheese_fungus.py \
    --hsi_path      /data/cheese/hyperspectral_cube.hdr \
    --csv_path      /home/a/openclaw/workspace/scripts/fungus_ref.csv \
    --output_dir    ./results \
    --thresh        0.13 \
    --mnf_components 18 \
    --min_mask_area 100 \
    --margin        3
```

Happy defect hunting! 🍄🧀  

---  

*End of README*  
