#!/usr/bin/env python3
"""Cellbin boundary/expression overlay template.

Template provenance:
Thor, DOI 10.1038/s41467-025-62593-1,
https://github.com/GuangyuWangLab2021/Thor, original files
`src/thor/plotting/boundary.py`, `spot_overlap.py`, and `spot_pie.py`.
BIDCell/CellSPA, DOI 10.1038/s41467-023-44560-w,
https://github.com/SydneyBioX/BIDCell and https://github.com/SydneyBioX/CellSPA,
original file `R/visualisation.R`.
Ascidian endostyle Stereo-seq, DOI 10.1126/sciadv.adi9035,
https://github.com/lskfs/ascidian-endostyle, original files
`01.stereo-seq_cellseg/cellsegment.py` and `gem_mask.py`.

Reusable success pattern: draw image-derived cell boundaries and expression/bin
coordinates in the same coordinate frame before interpreting cell-level labels.
Keep geometry QC separate from downstream clustering/DEG claims.
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


def require_polygon_columns(df: pd.DataFrame, id_col: str, x_col: str, y_col: str) -> None:
    missing = {id_col, x_col, y_col}.difference(df.columns)
    if missing:
        raise SystemExit(f"Boundary table is missing columns: {', '.join(sorted(missing))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundaries", required=True, help="CSV/TSV with one row per polygon vertex.")
    parser.add_argument("--expression", default="", help="Optional CSV/TSV with expression/bin/cell coordinates.")
    parser.add_argument("--boundary-id-col", default="cell_id")
    parser.add_argument("--boundary-x-col", default="x")
    parser.add_argument("--boundary-y-col", default="y")
    parser.add_argument("--expr-x-col", default="x")
    parser.add_argument("--expr-y-col", default="y")
    parser.add_argument("--expr-color-col", default="")
    parser.add_argument("--out", required=True)
    parser.add_argument("--point-size", type=float, default=0.6)
    parser.add_argument("--max-boundaries", type=int, default=5000)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--width", type=float, default=6.0)
    parser.add_argument("--height", type=float, default=5.6)
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    boundaries = read_table(args.boundaries)
    require_polygon_columns(boundaries, args.boundary_id_col, args.boundary_x_col, args.boundary_y_col)
    keep_ids = boundaries[args.boundary_id_col].drop_duplicates().head(args.max_boundaries)
    boundaries = boundaries[boundaries[args.boundary_id_col].isin(keep_ids)]

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    if args.expression:
        expr = read_table(args.expression)
        missing = {args.expr_x_col, args.expr_y_col}.difference(expr.columns)
        if missing:
            raise SystemExit(f"Expression table is missing columns: {', '.join(sorted(missing))}")
        if args.expr_color_col and args.expr_color_col in expr.columns:
            values = pd.to_numeric(expr[args.expr_color_col], errors="coerce")
            order = np.argsort(values.fillna(values.min()).to_numpy())
            sca = ax.scatter(
                expr[args.expr_x_col].to_numpy()[order],
                expr[args.expr_y_col].to_numpy()[order],
                c=values.to_numpy()[order],
                s=args.point_size,
                cmap="viridis",
                linewidths=0,
                alpha=0.75,
            )
            cbar = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
            cbar.ax.tick_params(labelsize=9)
        else:
            ax.scatter(
                expr[args.expr_x_col],
                expr[args.expr_y_col],
                s=args.point_size,
                c="#BAB0AC",
                linewidths=0,
                alpha=0.45,
                label="expression coordinates",
            )

    for _, sub in boundaries.groupby(args.boundary_id_col, sort=False):
        if len(sub) < 2:
            continue
        x = sub[args.boundary_x_col].to_numpy()
        y = sub[args.boundary_y_col].to_numpy()
        ax.plot(x, y, color="#E15759", linewidth=0.25, alpha=0.55)

    ax.set_aspect("equal")
    if args.invert_y:
        ax.invert_yaxis()
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.tick_params(labelsize=9)
    if not args.expr_color_col:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9, title_fontsize=10)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote cellbin boundary overlay to {out}")


if __name__ == "__main__":
    main()
