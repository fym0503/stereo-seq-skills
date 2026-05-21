#!/usr/bin/env python3
"""Paper-quality marker spatial grid for Stereo-seq h5ad objects.

Template provenance:
P09 AD prefrontal cortex, DOI 10.1038/s41467-024-54715-y,
Zenodo 10.5281/zenodo.14048103, original files `2_layer_annotation.ipynb`
and `2_Cell_bin_ref.ipynb`.
Human endometrium PCOS atlas, DOI 10.1038/s41591-025-03592-z,
https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R, original file
`Spatial RNA-seq analysis/04.1_Marker_Visualisation_Full_Slide.R`.
Avian optic tectum atlas, DOI 10.1016/j.isci.2024.109009,
https://github.com/Coleliao/Spatial_OT, original files
`spatial_analysis/fig1BC_spatial_visualization_mannaul_annotation.R` and
`spatial_analysis/fig2D_vlnplot_selected_deg.R`.

Reusable success pattern: show marker expression in equal-aspect tissue maps,
draw high-expression points on top, keep colorbars compact, use Arial-compatible
PDF output, and use one panel per gene to avoid legend overlap.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def require_scanpy():
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy`. Install or activate scanpy before "
            "running marker spatial grid plotting. Blocked figure: marker grid."
        ) from exc
    return sc


def parse_genes(value: str) -> list[str]:
    genes = [x.strip() for x in value.split(",") if x.strip()]
    if not genes:
        raise SystemExit("Provide at least one gene with --genes.")
    return genes


def gene_values(adata, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--genes", required=True, help="Comma-separated marker genes.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--ncols", type=int, default=4)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--clip-quantile", type=float, default=0.99)
    parser.add_argument("--width-per-panel", type=float, default=3.0)
    parser.add_argument("--height-per-panel", type=float, default=3.0)
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sc = require_scanpy()
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    genes = parse_genes(args.genes)
    missing = [gene for gene in genes if gene not in adata.var_names]
    if missing:
        print("Skip missing genes: " + ", ".join(missing))
    genes = [gene for gene in genes if gene in adata.var_names]
    if not genes:
        raise SystemExit("None of the requested genes are present in the AnnData object.")

    xy = np.asarray(adata.obsm[args.spatial_key])
    ncols = max(1, args.ncols)
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(args.width_per_panel * ncols, args.height_per_panel * nrows),
        squeeze=False,
    )
    for ax, gene in zip(axes.ravel(), genes):
        values = gene_values(adata, gene)
        vmax = float(np.nanquantile(values, args.clip_quantile))
        values = np.clip(values, 0, vmax if vmax > 0 else np.nanmax(values))
        order = np.argsort(values)
        sca = ax.scatter(
            xy[order, 0],
            xy[order, 1],
            c=values[order],
            s=args.point_size,
            cmap="magma",
            linewidths=0,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        if args.invert_y:
            ax.invert_yaxis()
        ax.set_title(gene, fontsize=11, pad=4)
        cbar = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.02)
        cbar.ax.tick_params(labelsize=8)
    for ax in axes.ravel()[len(genes) :]:
        ax.axis("off")
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote marker spatial grid to {out}")


if __name__ == "__main__":
    main()
