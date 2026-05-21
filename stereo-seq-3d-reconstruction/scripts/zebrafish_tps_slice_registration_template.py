#!/usr/bin/env python3
"""Register one moving Stereo-seq section to a reference mask with TPS.

Template provenance:
Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the
regenerating zebrafish heart
Paper DOI: 10.1038/s41467-025-59070-0
Original code: BGI-Qingdao/ZebrafishHeartRegeneration_project,
  10. Slice registrate to 3D model/CoornidateMatch_*.ipynb
Reusable success pattern: use tissue masks/anchors to transform slice
coordinates, then export transformed coordinates and an equal-aspect QC plot.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
try:
    from scipy.interpolate import Rbf
    from skimage import io, measure
except Exception as exc:  # pragma: no cover - environment guard
    raise SystemExit(
        "Missing Python packages `scipy` and/or `scikit-image`. Install them "
        "before running this template. Blocked step: TPS slice registration."
    ) from exc


def load_mask(path: str) -> np.ndarray:
    image = io.imread(path)
    if image.ndim == 3:
        image = image[..., :3].mean(axis=2)
    return image > np.nanmean(image)


def largest_contour(mask: np.ndarray) -> np.ndarray:
    contours = measure.find_contours(mask.astype(float), 0.5)
    if not contours:
        raise SystemExit("No tissue contour found in mask.")
    contour = max(contours, key=len)
    return np.column_stack([contour[:, 1], contour[:, 0]])


def sample_contour(contour: np.ndarray, n: int) -> np.ndarray:
    closed = np.vstack([contour, contour[0]])
    step = np.sqrt(((closed[1:] - closed[:-1]) ** 2).sum(axis=1))
    dist = np.concatenate([[0], np.cumsum(step)])
    targets = np.linspace(0, dist[-1], n, endpoint=False)
    x = np.interp(targets, dist, closed[:, 0])
    y = np.interp(targets, dist, closed[:, 1])
    return np.column_stack([x, y])


def fit_tps(moving: np.ndarray, reference: np.ndarray, smooth: float):
    fx = Rbf(moving[:, 0], moving[:, 1], reference[:, 0], function="thin_plate", smooth=smooth)
    fy = Rbf(moving[:, 0], moving[:, 1], reference[:, 1], function="thin_plate", smooth=smooth)
    return lambda points: np.column_stack([fx(points[:, 0], points[:, 1]), fy(points[:, 0], points[:, 1])])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-mask", required=True, help="Reference tissue mask image.")
    parser.add_argument("--moving-mask", required=True, help="Moving tissue mask image.")
    parser.add_argument("--moving-points", default="", help="CSV/TSV with moving x/y coordinates to transform.")
    parser.add_argument("--manual-anchors", default="", help="CSV/TSV with ref_x,ref_y,moving_x,moving_y columns.")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--n-anchors", type=int, default=80, help="Contour anchors when manual anchors are absent.")
    parser.add_argument("--smooth", type=float, default=0.0)
    parser.add_argument("--out-prefix", default="slice_registration")
    args = parser.parse_args()

    plt.rcParams["font.family"] = "Arial"
    out_prefix = Path(args.out_prefix)
    ref_mask = load_mask(args.ref_mask)
    mov_mask = load_mask(args.moving_mask)
    ref_contour = largest_contour(ref_mask)
    mov_contour = largest_contour(mov_mask)

    if args.manual_anchors:
        anchors = pd.read_csv(args.manual_anchors, sep=None, engine="python")
        needed = ["ref_x", "ref_y", "moving_x", "moving_y"]
        missing = [col for col in needed if col not in anchors.columns]
        if missing:
            raise SystemExit(f"Manual anchor file is missing columns: {missing}")
        ref_anchor = anchors[["ref_x", "ref_y"]].to_numpy(float)
        mov_anchor = anchors[["moving_x", "moving_y"]].to_numpy(float)
    else:
        ref_anchor = sample_contour(ref_contour, args.n_anchors)
        mov_anchor = sample_contour(mov_contour, args.n_anchors)
        anchors = pd.DataFrame(
            {
                "ref_x": ref_anchor[:, 0],
                "ref_y": ref_anchor[:, 1],
                "moving_x": mov_anchor[:, 0],
                "moving_y": mov_anchor[:, 1],
            }
        )

    transform = fit_tps(mov_anchor, ref_anchor, args.smooth)
    anchors.to_csv(out_prefix.with_suffix(".anchors.csv"), index=False)
    transformed_contour = transform(mov_contour)

    if args.moving_points:
        points = pd.read_csv(args.moving_points, sep=None, engine="python")
        for col in [args.x_col, args.y_col]:
            if col not in points.columns:
                raise SystemExit(f"Moving points file is missing column: {col}")
        xy = points[[args.x_col, args.y_col]].to_numpy(float)
        mapped = transform(xy)
        points["registered_x"] = mapped[:, 0]
        points["registered_y"] = mapped[:, 1]
        points.to_csv(out_prefix.with_suffix(".registered_points.csv"), index=False)

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.plot(ref_contour[:, 0], ref_contour[:, 1], color="#222222", lw=1.1, label="reference contour")
    ax.plot(transformed_contour[:, 0], transformed_contour[:, 1], color="#D55E00", lw=1.0, label="registered moving contour")
    ax.scatter(ref_anchor[:, 0], ref_anchor[:, 1], s=10, color="#0072B2", label="reference anchors", zorder=3)
    ax.scatter(transform(mov_anchor)[:, 0], transform(mov_anchor)[:, 1], s=8, color="#CC79A7", label="mapped moving anchors", zorder=3)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".qc.pdf"), dpi=300)
    plt.close(fig)
    print(f"Wrote registration outputs with prefix {out_prefix}")


if __name__ == "__main__":
    main()
