# 📊 TESTING README – Cheese Fungus Segmentation Integration Test

**Location:** `workspace/scripts/TESTING_README.md`  
**Purpose:** Guide for running the end‑to‑end integration test that validates the synthetic hyperspectral cube, processing script, and output artefacts.

---  

## 📂 Directory Layout (relevant files)

| File | Description |
|------|-------------|
| `generate_test_data.py` | Creates `test_cube.npy` (synthetic HSI) and `test_fungus_ref.csv` (reference spectra). |
| `segment_cheese_fungus_integration.py` | Processing script that reads the synthetic data, runs MNF → NMF → threshold → cleaning, and writes mask, overlay, and summary CSV. |
| `test_integration.py` | Pytest‑style test that invokes the processing script and asserts that all three output artefacts are created and non‑empty. |
| `requirements.txt` | Minimal list of Python packages required (`numpy`, `scipy`, `rasterio`, `matplotlib`, `tqdm`, `scikit-learn`). |
| `TESTING_README.md` | This document – step‑by‑step instructions for executing the test suite. |

---  

## 🛠️ Prerequisites – Install Required Packages

The integration test relies on a few third‑party libraries that may not be present in a stock Python installation.

### Option A – System‑wide installation (quick start)

```bash
# From the workspace root (where the `scripts` folder lives)
python -m pip install -r scripts/requirements.txt --break-system-packages
```

### Option B – Isolated virtual environment (recommended for reproducibility)

```bash
# 1️⃣ Create a clean virtual environment
python -m venv .venv          # creates .venv/ in the workspace root

# 2️⃣ Activate it
source .venv/bin/activate      # Bash/Zsh
# Windows PowerShell: .venv\Scripts\Activate.ps1

# 3️⃣ Install the dependencies
python -m pip install -r scripts/requirements.txt
```

> **Tip:** Activating the environment ensures that `python` now refers to the isolated interpreter with all required packages available.

---  

## 📦 Step‑by‑Step Execution Workflow

Below is the exact command sequence you can copy‑paste into your terminal.  
All commands are relative to the **workspace root** (the folder that contains the `scripts/` directory).

```bash
# --------------------------------------------------------------
# 1️⃣  Generate the synthetic test data (one‑time only)
# --------------------------------------------------------------
cd workspace/scripts
python generate_test_data.py
# <-- Outputs: test_cube.npy  and  test_fungus_ref.csv

# --------------------------------------------------------------
# 2️⃣  Install the required Python packages (if not already done)
# --------------------------------------------------------------
python -m pip install -r requirements.txt --break-system-packages
# --------------------------------------------------------------
# 3️⃣  Run the processing script on the synthetic data
# --------------------------------------------------------------
python segment_cheese_fungus_integration.py \
    --hsi_path test_cube.npy \
    --csv_path test_fungus_ref.csv \
    --output_dir integration_test_output \
    --thresh 0.12 \
    --mnf_components 5 \
    --min_mask_area 5 \
    --margin 2
# <-- Creates: integration_test_output/test_cube_mask.tif
#                integration_test_output/test_cube_overlay.png
#                integration_test_output/test_cube_summary.csv

# --------------------------------------------------------------
# 4️⃣  Validate the output with the pytest integration test
# --------------------------------------------------------------
python -m pytest test_integration.py -v
# <-- Expected output:  test_integration PASSED
```

---  

## ✅ What the Test Verifies

| Assertion | Reason |
|-----------|--------|
| **Script exits with code 0** | No uncaught exceptions during processing. |
| **All three output files exist** | Mask (GeoTIFF), overlay (PNG), and summary (CSV) are produced. |
| **Files are non‑empty** | Guarantees that data was actually written (size > 0 bytes). |
| **CSV contains numeric values ≥ 0** | Sanity‑check that the summary statistics are meaningful. |

If any assertion fails, the test will raise an `AssertionError` and print a helpful message with the offending file/column.

---  

## 🛠️ Debugging & Common Issues

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `ModuleNotFoundError: No module named 'rasterio'` | `rasterio` not installed | `python -m pip install rasterio --break-system-packages` |
| `ImportError: cannot import name 'FactorAnalysis'` | `scikit‑learn` missing | `python -m pip install scikit-learn --break-system-packages` |
| Empty output files | Threshold too high or MNF components insufficient | Reduce `--thresh` (e.g., `0.05`) or increase `--mnf_components`. |
| Test fails with “File not found” | Running the pytest command from the wrong directory | `cd workspace/scripts` before invoking `python -m pytest`. |
| Permission denied when installing packages | System‑wide pip not allowed | Use a virtual environment (highly recommended) or add `--break-system-packages`. |

---  

## 🧪 Extending the Test Suite

You can tweak any of the command‑line arguments to explore how they affect the results:

| Argument | What it Controls | Typical Experiments |
|----------|------------------|---------------------|
| `--thresh` | Abundance cut‑off for binary mask | Test `0.05`, `0.10`, `0.15` |
| `--mnf_components` | Number of MNF components kept | Try `3`, `10`, `20` |
| `--min_mask_area` | Minimum object size (pixels) after cleaning | Test `2`, `10`, `20` |
| `--margin` | Pixel margin added to the mask for the overlay | Try `0`, `5`, `10` |
| `--csv_path` | Reference CSV file | Replace with a real fungal spectrum if available |

After changing an argument, simply re‑run steps **3️⃣** and **4️⃣** to see the effect.

---  

## 📚 Citations & Acknowledgements

- **CoreSpecViewer** – Russell Rogers, Geological Survey Ireland.  
- **Spectral Python (SPy)** – https://github.com/spectralpython/spectral  
- **OpenClaw** – Open‑source AI assistant framework (https://openclaw.ai)  

If you use this test suite in a publication, please cite the above projects where appropriate.

---  

## 🎉 You’re Ready to Test!

Run the commands in the **“Step‑by‑Step Execution Workflow”** section, watch the progress messages, and verify that the final pytest output reads:

```
✅ test_integration PASSED
✅ All integration-test checks passed!
```

That means the synthetic data was processed successfully, the three required artefacts were generated, and the automated validation confirmed everything looks good.

Happy testing! 🚀  

---  

*End of TESTING_README.md*  
