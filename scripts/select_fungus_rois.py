#!/usr/bin/env python3
"""
Interactive ROI selector for fungus regions in HSI cubes.

Usage (from workspace root or scripts folder, with venv active):

    python scripts/select_fungus_rois.py \
        --hsi_path "D:\PycharmProjects\cube_26_02_12_37_27\cube_26_02_12_37_27.hdr" \
        --output_mask "D:\PycharmProjects\Workspace\roi_fungus_mask.npy"

Controls:
    - Left‑drag   : draw rectangle ROI
    - Right‑drag  : draw ellipse (circle/oval) ROI
    - 'z' key     : clear all ROIs
    - 's' key     : save mask and exit
    - 'q' / ESC   : exit without saving
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector, EllipseSelector


def load_hsi_for_display(hsi_path: str):
    """
    Load HSI cube and derive a simple RGB preview for interaction.
    Supports ENVI .hdr/.dat via spectral.io.envi and .npy stacks.
    """
    lower = hsi_path.lower()
    if lower.endswith((".hdr", ".dat")):
        import spectral.io.envi as envi

        img = envi.open(hsi_path)
        cube = img.load().astype(np.float32)  # (rows, cols, bands)
    else:
        cube = np.load(hsi_path).astype(np.float32)

    h, w, b = cube.shape

    # Choose three bands spread across the spectrum for a rough RGB
    b1 = 0
    b2 = b // 2
    b3 = max(b - 1, 0)
    rgb = np.stack([cube[:, :, b1], cube[:, :, b2], cube[:, :, b3]], axis=-1)
    # Normalize to [0, 1] for display
    rgb_min = np.percentile(rgb, 2)
    rgb_max = np.percentile(rgb, 98)
    rgb = np.clip((rgb - rgb_min) / (rgb_max - rgb_min + 1e-6), 0, 1)
    return cube, rgb


def main():
    parser = argparse.ArgumentParser(
        description="Interactively select rectangular and circular ROIs on an HSI cube."
    )
    parser.add_argument(
        "--hsi_path",
        required=True,
        help="Path to HSI cube (.hdr/.dat ENVI or .npy stack).",
    )
    parser.add_argument(
        "--output_mask",
        required=True,
        help="Output .npy file for boolean ROI mask (shape: rows x cols).",
    )
    args = parser.parse_args()

    cube, rgb = load_hsi_for_display(args.hsi_path)
    h, w, _ = rgb.shape
    roi_mask = np.zeros((h, w), dtype=bool)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title("ROI selection: left‑drag=rect, right‑drag=ellipse, 's'=save, 'z'=clear, 'q'=quit")
    img_artist = ax.imshow(rgb)
    overlay_artist = ax.imshow(
        np.zeros((h, w)), cmap="Reds", alpha=0.0
    )  # initially invisible overlay

    def update_overlay():
        # Update overlay data in red where roi_mask is True
        if roi_mask.any():
            roi_vis = np.ma.masked_where(~roi_mask, roi_mask)
            overlay_artist.set_data(roi_vis)
            overlay_artist.set_alpha(0.4)
        else:
            overlay_artist.set_alpha(0.0)
        fig.canvas.draw_idle()

    def onselect_rect(eclick, erelease):
        if eclick.xdata is None or erelease.xdata is None:
            return
        x0, y0 = int(round(eclick.xdata)), int(round(eclick.ydata))
        x1, y1 = int(round(erelease.xdata)), int(round(erelease.ydata))
        x_min, x_max = sorted((x0, x1))
        y_min, y_max = sorted((y0, y1))
        x_min = max(x_min, 0)
        y_min = max(y_min, 0)
        x_max = min(x_max, w - 1)
        y_max = min(y_max, h - 1)
        roi_mask[y_min : y_max + 1, x_min : x_max + 1] = True
        update_overlay()

    def onselect_ellipse(eclick, erelease):
        if eclick.xdata is None or erelease.xdata is None:
            return
        x0, y0 = eclick.xdata, eclick.ydata
        x1, y1 = erelease.xdata, erelease.ydata
        cx = 0.5 * (x0 + x1)
        cy = 0.5 * (y0 + y1)
        rx = abs(x1 - x0) / 2.0
        ry = abs(y1 - y0) / 2.0
        if rx < 1 or ry < 1:
            return
        yy, xx = np.ogrid[:h, :w]
        mask_ellipse = (((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2) <= 1.0
        roi_mask[mask_ellipse] = True
        update_overlay()

    rect_selector = RectangleSelector(
        ax,
        onselect_rect,
        useblit=True,
        button=[1],  # left mouse
        minspanx=2,
        minspany=2,
        spancoords="pixels",
    )

    ellipse_selector = EllipseSelector(
        ax,
        onselect_ellipse,
        useblit=True,
        button=[3],  # right mouse
        minspanx=2,
        minspany=2,
        spancoords="pixels",
    )

    saved = {"done": False}

    def on_key(event):
        if event.key == "z":
            roi_mask[:] = False
            update_overlay()
        elif event.key == "s":
            os.makedirs(os.path.dirname(args.output_mask) or ".", exist_ok=True)
            np.save(args.output_mask, roi_mask.astype(bool))
            print(f"ROI mask saved to {args.output_mask}")
            saved["done"] = True
            plt.close(fig)
        elif event.key in ("q", "escape"):
            print("Exiting without saving.")
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)

    update_overlay()
    plt.show()

    if not saved["done"]:
        # If user closed the window without pressing 's', do not save.
        pass


if __name__ == "__main__":
    main()

