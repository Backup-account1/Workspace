## Fungus ROI Selection Tool (`select_fungus_rois.py`)

This tool lets you **manually select fungus regions** in a hyperspectral image (HSI) by drawing rectangles and circles/ellipses, and saves the result as a boolean ROI mask that can be used by the segmentation script.

---

### 1. Purpose

- **Visually mark** regions where fungus is present (or where you want to focus analysis).
- Generate a **2‑D boolean mask** (`rows x cols`) as a `.npy` file.
- Use that mask to **restrict segmentation** so only pixels inside your ROIs are kept.

---

### 2. Requirements

- Python running in your existing workspace venv:
  - `matplotlib`
  - `numpy`
  - `spectral` (for ENVI `.hdr/.dat` support)

These are already installed in `venv` from the earlier setup.

---

### 3. Running the ROI Tool

#### 3.1 Generic example

From the workspace root (`D:\PycharmProjects\Workspace`):

```bash
venv\Scripts\python scripts\select_fungus_rois.py ^
  --hsi_path "D:\PycharmProjects\cube_26_02_12_37_27\cube_26_02_12_37_27.hdr" ^
  --output_mask "D:\PycharmProjects\Workspace\roi_fungus_mask.npy"
```

#### 3.2 Using an HSI from the `capture` folder

If your ENVI cube lives in the workspace `capture` directory, run (still from the workspace root):

```bash
venv\Scripts\python scripts\select_fungus_rois.py ^
  --hsi_path "D:\PycharmProjects\Workspace\capture\2025-10-31_09-03-44_white_circ.hdr" ^
  --output_mask "D:\PycharmProjects\Workspace\capture\2025-10-31_09-03-44_white_circ_roi_mask.npy"
```

You can swap in any other `.hdr` in `capture` and change the mask filename as needed.

**Arguments**

- `--hsi_path`  
  Path to the HSI cube:
  - ENVI: `.hdr` (with its paired `.dat`)  
  - or `.npy` stack with shape `(rows, cols, bands)`

- `--output_mask`  
  Where to save the ROI mask as a `.npy` file.  
  The saved array will have shape `(rows, cols)` and `dtype=bool`.

---

### 4. Mouse & Keyboard Controls

In the interactive window:

- **Left mouse drag**: draw a **rectangular** ROI
- **Right mouse drag**: draw a **circular/elliptical** ROI
- **Mouse wheel**: zoom in/out around the cursor
- **z**: clear all ROIs (reset mask to all‑False)
- **s**: save current mask to `--output_mask` and exit
- **q** or **Esc**: quit **without** saving

Multiple rectangles/circles are **combined** (logical OR) into a single mask.

---

### 5. How the Display Works

- The script loads the full cube and builds a simple **RGB preview** by taking:
  - First band
  - Middle band
  - Last band
- It stretches intensities between the 2‑nd and 98‑th percentiles so the image is visible.
- The red overlay shows the current ROI mask.

This is only for visualization – the full cube is still available for later processing.

---

### 6. Using the ROI Mask in Segmentation

Once you’ve saved `roi_fungus_mask.npy`, run the segmentation script with `--roi_mask`:

```bash
cd D:\PycharmProjects\Workspace\scripts

..\venv\Scripts\python segment_cheese_fungus_integration.py ^
  --hsi_path "D:\PycharmProjects\cube_26_02_12_37_27\cube_26_02_12_37_27.hdr" ^
  --csv_path "D:\PycharmProjects\Workspace\test_fungus_ref.csv" ^
  --output_dir "D:\PycharmProjects\Workspace\real_cube_26_02_12_37_27_output_roi" ^
  --thresh 0.12 ^
  --mnf_components 5 ^
  --min_mask_area 5 ^
  --margin 2 ^
  --roi_mask "D:\PycharmProjects\Workspace\roi_fungus_mask.npy"
```

Behavior:

- The pipeline runs as usual (MNF → NMF → threshold → cleaning).
- At the end, the cleaned mask is **AND‑ed** with your ROI mask:
  - Pixels outside the ROI become 0.
  - Pixels inside the ROI keep their predicted fungus/non‑fungus label.

Resulting files (in `--output_dir`):

- `*_mask.tif` – binary mask limited to your ROIs
- `*_overlay.png` – RGB preview with both segmentation and ROI applied
- `*_summary.csv` – statistics for pixels inside your ROIs only

---

### 7. Tips

- Draw ROIs **only over regions you want to analyse** (e.g., cheese surface with visible fungus).
- Use multiple small rectangles/circles to avoid including irrelevant background.
- If the mask shape mismatch error appears, make sure:
  - `--hsi_path` in the ROI tool and segmentation command point to the **same cube**.
  - You didn’t crop or resample the cube between the two steps.

