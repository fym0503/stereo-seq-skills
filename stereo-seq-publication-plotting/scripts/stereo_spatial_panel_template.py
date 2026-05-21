#!/usr/bin/env python3
"""Paper-quality Stereo-seq spatial map template.

Template provenance:
P09 AD prefrontal cortex, DOI 10.1038/s41467-024-54715-y,
Zenodo 10.5281/zenodo.14048103, original files `2_layer_annotation.ipynb`
and `2_Cell_bin_ref.ipynb`.
Human endometrium PCOS atlas, DOI 10.1038/s41591-025-03592-z,
https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R, original file
`Spatial RNA-seq analysis/00.1_Stereopy_binning.py`.
GF/SPF cecum atlas, DOI 10.1016/j.isci.2024.108941,
https://github.com/1014723815/GF_SPF_cecum, original file
`Spatial transcriptome.R`.

Reusable success pattern: equal-aspect tissue maps, stable categorical colors,
readable legends outside the data region, Arial-compatible PDF output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_PALETTE = [
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#B07AA1",
    "#76B7B2",
    "#EDC948",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#1F77B4",
    "#D62728",
    "#2CA02C",
    "#9467BD",
    "#8C564B",
    "#17BECF",
]

P09_LAYER_PALETTE = {
    "L1": "#FF7F00",
    "L2/3": "#984EA3",
    "L4": "#4DAF4A",
    "L5": "#377EB8",
    "L6": "#E41A1C",
    "WM": "#A65628",
}


def require_scanpy():
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy`. Install or activate an environment "
            "with scanpy before reading .h5ad input. Blocked step: spatial plotting."
        ) from exc
    return sc


def read_input(path: str, x_col: str, y_col: str, color_col: str, spatial_key: str) -> pd.DataFrame:
    path_obj = Path(path)
    if path_obj.suffix.lower() == ".h5ad":
        sc = require_scanpy()
        adata = sc.read_h5ad(path)
        if spatial_key not in adata.obsm:
            raise SystemExit(f"AnnData is missing obsm['{spatial_key}'].")
        if color_col not in adata.obs and color_col not in adata.var_names:
            raise SystemExit(f"`{color_col}` is not in adata.obs or adata.var_names.")
        xy = adata.obsm[spatial_key]
        df = pd.DataFrame({x_col: xy[:, 0], y_col: xy[:, 1]}, index=adata.obs_names)
        if color_col in adata.obs:
            df[color_col] = adata.obs[color_col].astype(str).to_numpy()
        else:
            values = adata[:, color_col].X
            if hasattr(values, "toarray"):
                values = values.toarray()
            df[color_col] = np.asarray(values).ravel()
        return df
    return pd.read_csv(path, sep=None, engine="python")


def build_palette(values: list[str], palette_path: str | None) -> dict[str, str]:
    if palette_path:
        tab = pd.read_csv(palette_path, sep=None, engine="python")
        if not {"label", "color"}.issubset(tab.columns):
            raise SystemExit("Palette table must contain columns: label, color")
        return dict(zip(tab["label"].astype(str), tab["color"].astype(str)))
    palette = dict(P09_LAYER_PALETTE)
    missing = [v for v in values if v not in palette]
    for i, value in enumerate(missing):
        palette[value] = BASE_PALETTE[i % len(BASE_PALETTE)]
    return palette


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CSV/TSV table or .h5ad")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--color-col", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--palette", default=None, help="Optional table with label,color columns")
    parser.add_argument("--out", required=True, help="Output PDF/PNG path")
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--alpha", type=float, default=0.95)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--continuous", action="store_true")
    parser.add_argument("--title", default="")
    parser.add_argument("--width", type=float, default=6.0)
    parser.add_argument("--height", type=float, default=5.4)
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    df = read_input(args.input, args.x_col, args.y_col, args.color_col, args.spatial_key)
    required = {args.x_col, args.y_col, args.color_col}
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"Input is missing columns: {', '.join(sorted(missing))}")
    df = df.dropna(subset=[args.x_col, args.y_col, args.color_col]).copy()

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    if args.continuous or pd.api.types.is_numeric_dtype(df[args.color_col]):
        values = pd.to_numeric(df[args.color_col], errors="coerce")
        order = np.argsort(values.to_numpy())
        sca = ax.scatter(
            df[args.x_col].to_numpy()[order],
            df[args.y_col].to_numpy()[order],
            c=values.to_numpy()[order],
            s=args.point_size,
            cmap="viridis",
            alpha=args.alpha,
            linewidths=0,
        )
        cbar = fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.02)
        cbar.ax.tick_params(labelsize=9)
    else:
        labels = sorted(df[args.color_col].astype(str).unique())
        palette = build_palette(labels, args.palette)
        for label in labels:
            sub = df[df[args.color_col].astype(str) == label]
            ax.scatter(
                sub[args.x_col],
                sub[args.y_col],
                s=args.point_size,
                color=palette[label],
                label=label,
                alpha=args.alpha,
                linewidths=0,
            )
        legend_cols = 1 if len(labels) <= 18 else 2
        leg = ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
            frameon=False,
            markerscale=max(2.5, 5 / max(args.point_size, 0.5)),
            fontsize=9,
            title=args.color_col,
            title_fontsize=10,
            ncol=legend_cols,
        )
        legend_handles = getattr(leg, "legend_handles", getattr(leg, "legendHandles", []))
        for handle in legend_handles:
            handle.set_alpha(1)

    ax.set_aspect("equal")
    ax.axis("off")
    if args.invert_y:
        ax.invert_yaxis()
    if args.title:
        ax.set_title(args.title, fontsize=11, pad=6)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote spatial panel to {out}")


if __name__ == "__main__":
    main()
