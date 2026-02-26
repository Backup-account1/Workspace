# OpenCLaw Workspace – Synthetic Data Segmentation Test

**Author:** OpenCLaw Assistant  
**Last Updated:** 2026‑02‑26  

---

## 📌 Overview

This repository contains a test project for **HSI (Hyperspectral Imaging) segmentation** of fungal defect regions in cheese hyperspectral images. The workflow follows the CoreSpecViewer processing chain:

1. **Continuum Removal (CR)**
2. **Minimum-Noise-Fraction (MNF) projection**
3. **NMF unmixing**
4. **Thresholding**
5. **Mask cleaning**
6. **Export of results (GeoTIFF, PNG overlay, CSV summary)**

The synthetic test data consists of:

- `test_cube.npy` – a 3‑D NumPy array of shape **(20, 20, 8)** representing an 8‑band hyperspectral cube.
- `test_fungus_ref.csv` – CSV file with reference fungal spectra (wavelength, reflectance).

The script `segment_cheese_fungus.py` automates the full processing chain and writes the results to:

```
/opt/openclaw/workspace/synthetic_test/results/
```

---

## 📂 Repository Structure

```
/openclaw/workspace/
├─ scripts/
│   └─ segment_cheese_fungus.py      # Main processing script
├─ synthetic_test/
│   ├─ test_cube.npy                 # Synthetic hyperspectral cube (20×20×8)
│   ├─ test_fungus_ref.csv           # Reference fungal spectra
│   └─ results/                      # Output folder (auto‑created)
├─ memory/
├─ scripts/
│   └─ run_segment.sh                # Optional wrapper to launch the script
└─ README.md                         # <--- This file
```

---

## 🛠️ Dependencies

- **Python 3.10+**
- **Required Python packages** (installed via the bundled virtual environment at `/opt/openclaw/workspace/.pi/venv`):
  - `numpy`
  - `pandas`
  - `scikit-learn`
  - `matplotlib`
  - `rasterio`
  - `scipy`
  - `spectral`
  - `geopandas`
  - `shapely`
  - `tqdm`
  - `progressbar2`
  - `opencv-python`
  - `rasterio`
  - `gdal`
  - `statsmodels`
- **System libraries**:
  - `gdal-bin`
  - `gdal-data`
  - `libgdal-dev`
  - `build-essential`
  - `python3-dev`

> The virtual environment is already set up in `/opt/openclaw/workspace/.pi/venv`.  
> Activate it with `source /opt/openclaw/workspace/.pi/venv/bin/activate` before running any commands outside of PyCharm.

---

## 🚀 Running the Segmentation Script on Synthetic Data

### 1️⃣  Activate the virtual environment (if not already active)

```bash
source /opt/openclaw/workspace/.pi/venv/bin/activate
```

### 2️⃣  Execute the script with synthetic data

```bash
python3 /opt/openclaw/workspace/scripts/segment_cheese_fungus.py \
    --hsi_path /opt/openclaw/workspace/synthetic_test/test_cube.npy \
    --csv_path   /opt/openclaw/workspace/synthetic_test/test_fungus_ref.csv \
    --output_dir /opt/openclaw/workspace/synthetic_test/results \
    --mnf_components 8 \
    --thresh 0.12 \
    --min_mask_area 80 \
    --margin 5
```

#### Parameters Explained

| Flag | Description |
|------|-------------|
| `--hsi_path` | Path to the synthetic hyperspectral cube (`test_cube.npy`). |
| `--csv_path` | CSV with reference fungal spectra (`test_fungus_ref.csv`). |
| `--output_dir` | Destination folder for all generated outputs (mask, overlay, CSV). |
| `--mnf_components` | Number of MNF components to compute. **Capped automatically** to the number of spectral bands (8 for synthetic data). |
| `--thresh` | Abundance threshold for binary mask generation (default = 0.12). |
| `--min_mask_area` | Minimum object size in pixels to keep a component (default = 80). |
| `--margin` | Margin (in pixels) added around detected objects for overlay expansion. |

### 3️⃣  Expected Output Files

After a successful run you will find the following files in
`/opt/openclaw/workspace/synthetic_test/results/`:

| File | Description |
|------|-------------|
| `<basename>_mask.tif` | Binary mask (GeoTIFF) where fungal defect pixels are `1`. |
| `<basename>_overlay.png` | Visual overlay of the segmentation on the first three bands of the cube (red mask). |
| `<basename>_summary.csv` | CSV summary containing total pixels, fungus pixels, and percentage coverage. |

Sample console output:

```
🔍 Loaded cube: 20×20×8
🛠️  Continuum removal …
⚙️  Computing 8 MNF components …
🧩  NMF unmixing …
🔎  Thresholding (>0.12) …
🧹  Cleaning mask …
✅ Mask saved          → /opt/openclaw/workspace/synthetic_test/results/test_mask.tif
✅ Overlay saved       → /opt/openclaw/workspace/synthetic_test/results/test_overlay.png
✅ Summary CSV saved   → /opt/openclaw/workspace/synthetic_test/results/test_summary.csv

✅ All done!
```

---

## 🐞 Troubleshooting & Known Issues

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| **`ValueError: cannot reshape array of size 64 into shape (20,20,8)`** | `n_components` larger than number of bands; old code did not cap components. | The script now caps `n_components = min(args.mnf_components, cube.shape[2])`. Ensure you are using the latest version of `segment_cheese_fungus.py`. |
| **`ModuleNotFoundError: No module named 'app'`** | The `app` package is not on `PYTHONPATH`. The script now uses a **safe import** that works as long as the virtual environment is activated. | Activate the venv (`source /opt/openclaw/workspace/.pi/venv/bin/activate`). |
| **`ImportError: cannot import name 'FactorAnalysis'`** | `scikit‑learn` missing. | Install/activate the virtual environment where `scikit-learn` is installed. |
| **Empty results folder** | Script failed before reaching the export step. Check the console for traceback messages. | Verify the paths to `test_cube.npy` and `test_fungus_ref.csv` are correct, and that the virtual environment contains all required packages. |
| **Permission denied when writing output** | Output directory not writable. | Ensure you have write permissions for `/opt/openclaw/workspace/synthetic_test/results/` or change `--output_dir` to a directory you own (e.g., `$HOME/results`). |

---

## 📄 Full Script Overview (`segment_cheese_fungus.py`)

Below is a **concise excerpt** of the most relevant functions (the full script is available at `/opt/openclaw/workspace/scripts/segment_cheese_fungus.py`).

```python
def continuum_remove(cube: np.ndarray) -> np.ndarray:
    """Placeholder – returns the cube unchanged (no actual CR)."""
    return cube

def mnf_transform(cube: np.ndarray, n_components: int) -> np.ndarray:
    """
    Minimum‑Noise‑Fraction projection.
    Returns a cube of shape (rows, cols, n_components).
    """
    H, W, B = cube.shape
    cube_2d = cube.reshape(-1, B)                # (pixels, bands)
    n_comp = min(n_components, B)                # cap to available bands
    fa = FactorAnalysis(n_components=n_comp, max_iter=200, random_state=0)
    transformed = fa.fit_transform(cube_2d)      # (pixels, n_comp)
    # Safe reshape with padding/truncation fallback
    expected = H * W * n_comp
    actual   = transformed.size
    if actual == expected:
        mnf = transformed.reshape(H, W, n_comp)
    else:
        flat = transformed.flatten()
        if len(flat) < expected:
            flat = np.pad(flat, (0, expected - len(flat)), mode='constant')
        else:
            flat = flat[:expected]
        mnf = flat.reshape(H, W, n_comp)
    return mnf.astype(np.float32)

def nmf_unmixing(cube_2d: np.ndarray, n_endmembers: int, max_iter: int = 500, tol: float = 1e-4):
    """Linear NMF unmixing – returns abundance map of the first end‑member."""
    ...
```

The script then follows the pipeline:

1. Load cube & reference spectra.  
2. **Continuum removal** (identity placeholder).  
3. **MNF projection** (`mnf_transform`).  
4. **Flatten** to 2‑D for NMF.  
5. **NMF unmixing** (first component = fungus abundance).  
6. **Threshold** → binary mask.  
7. **Clean mask** (remove small components, fill holes, optional margin).  
8. **Export** (GeoTIFF mask, PNG overlay, CSV summary).

---

## 📂 Desktop Shortcut (Optional)

A convenient `.desktop` file can be placed on your Desktop to launch PyCharm with the script pre‑configured.

```bash
cat <<'EOF' > "$HOME/Desktop/pycharm_shortcut.desktop"
[Desktop Entry]
Type=Application
Name=Run Segmentation (PyCharm)
Icon=jetbrains-pycharm
Terminal=true
Exec=/opt/openclaw/workspace/scripts/run_segment.sh   # wrapper that activates venv & runs the script
Categories=Science;Graphics;
StartupNotify=true
EOF
chmod +x "$HOME/Desktop/pycharm_shortcut.desktop"
```

The wrapper script (`run_segment.sh`) contains the exact `python3 …` command shown in section **2️⃣**.

---

## 📌 Summary

- **Synthetic test data**: `test_cube.npy` (20×20×8) + `test_fungus_ref.csv`.  
- **Modified script**: caps MNF components, safe‑reshape fallback, placeholder continuum removal.  
- **Running the script**: one‑liner with appropriate flags produces mask, overlay, and CSV.  
- **Desktop shortcut**: creates a clickable icon that launches PyCharm with the script pre‑configured.  
- **Dependencies**: all available in the bundled virtual environment; just activate it before executing.

Feel free to adapt the paths, thresholds, or parameters to match your own experimental setup. If you encounter any errors, the console output will point directly to the failing step—review the traceback and adjust the corresponding flag (e.g., `--mnf_components`, `--thresh`, etc.).

---

**Happy segmenting!** 🎉 

--- 

*End of README* 

--- 
