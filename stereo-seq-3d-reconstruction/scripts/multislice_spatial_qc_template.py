#!/usr/bin/env python3
"""Multi-slice spatial coordinate QC template for Stereo-seq 3D work.

Template provenance:
SLAT, DOI 10.1038/s41467-023-43105-5,
https://github.com/gao-lab/SLAT, original files `benchmark/analysis/plot_slices.ipynb`,
`benchmark/analysis/plot_keypoints.ipynb`, and `benchmark/analysis/3d_analysis.ipynb`.
SPACEL, DOI 10.1038/s41467-023-43220-3,
https://github.com/QuKunLab/SPACEL, original files
`SPACEL/Scube/plot.py` and `docs/tutorials/Stereo-seq_Scube.ipynb`.
Zebrafish heart regeneration, DOI 10.1038/s41467-025-59070-0,
https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project, original
`10. Slice registrate to 3D model/CoornidateMatch_*.ipynb`.

Reusable success pattern: before 3D biological interpretation, show every slice
in a common coordinate style, export coordinate bounds/centroids, and plot a
3D-ready stacked scatter with equal 2D slice aspect.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coords", required=True, help="CSV/TSV with x, y, and slice columns.")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--slice-col", default="slice")
    parser.add_argument("--color-col", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="multislice_qc")
    parser.add_argument("--point-size", type=float, default=0.5)
    parser.add_argument("--slice-spacing", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    df = read_table(args.coords)
    missing = {args.x_col, args.y_col, args.slice_col}.difference(df.columns)
    if missing:
        raise SystemExit(f"Coordinate table is missing columns: {', '.join(sorted(missing))}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df[args.slice_col] = df[args.slice_col].astype(str)
    slices = list(dict.fromkeys(df[args.slice_col].tolist()))

    summary = df.groupby(args.slice_col).agg(
        n_units=(args.x_col, "size"),
        x_min=(args.x_col, "min"),
        x_max=(args.x_col, "max"),
        y_min=(args.y_col, "min"),
        y_max=(args.y_col, "max"),
        x_centroid=(args.x_col, "mean"),
        y_centroid=(args.y_col, "mean"),
    )
    summary.to_csv(outdir / f"{args.prefix}_slice_bounds.tsv", sep="\t")

    ncols = min(4, max(1, len(slices)))
    nrows = int(np.ceil(len(slices) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 3.0 * nrows), squeeze=False)
    palette = dict(zip(slices, plt.cm.tab20(np.linspace(0, 1, max(len(slices), 1)))))
    for ax, slice_name in zip(axes.ravel(), slices):
        sub = df[df[args.slice_col] == slice_name]
        if args.color_col and args.color_col in sub.columns:
            values = pd.to_numeric(sub[args.color_col], errors="coerce")
            ax.scatter(sub[args.x_col], sub[args.y_col], c=values, s=args.point_size, cmap="viridis", linewidths=0)
        else:
            ax.scatter(sub[args.x_col], sub[args.y_col], c=[palette[slice_name]], s=args.point_size, linewidths=0)
        ax.set_aspect("equal")
        ax.axis("off")
        if args.invert_y:
            ax.invert_yaxis()
        ax.set_title(slice_name, fontsize=10)
    for ax in axes.ravel()[len(slices) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / f"{args.prefix}_slice_panels.pdf", dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(6.2, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    for z, slice_name in enumerate(slices):
        sub = df[df[args.slice_col] == slice_name]
        ax.scatter(
            sub[args.x_col],
            sub[args.y_col],
            np.full(len(sub), z * args.slice_spacing),
            s=args.point_size,
            color=palette[slice_name],
            alpha=0.75,
            linewidths=0,
            label=slice_name,
        )
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.set_zlabel("slice", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9, title_fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f"{args.prefix}_stacked_3d_preview.pdf", dpi=300)
    plt.close(fig)
    print(f"Wrote multi-slice QC outputs to {outdir}")


if __name__ == "__main__":
    main()
