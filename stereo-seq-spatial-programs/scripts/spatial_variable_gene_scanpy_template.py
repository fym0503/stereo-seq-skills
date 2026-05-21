#!/usr/bin/env python3
"""Spatially variable gene/program template for Stereo-seq h5ad objects.

Template provenance:
SVGbench, DOI 10.1186/s13059-023-03145-y,
https://github.com/PYangLab/SVGbench, original files under
`02_run methods/run_spatialde/`.
PROST, DOI 10.1038/s41467-024-44835-w,
https://github.com/Tang-Lab-super/PROST, original file `test/Stereo-seq.ipynb`
and plotting helpers in `PROST/plot.py`.
SpaSEG, DOI 10.1186/s13059-025-03697-1,
https://github.com/y-bai/SpaSEG, original Stereo-seq marker/SVG notebooks and
`downstream/plotting/_svgplot.py`.

Reusable success pattern: identify spatially patterned genes with explicit
coordinates, export a ranked table, and immediately plot top genes in tissue
coordinates as evidence for domain/program interpretation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def require_packages():
    try:
        import scanpy as sc
        import squidpy as sq
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy` and/or `squidpy`. Install them before "
            "running spatial autocorrelation. Blocked step: spatially variable gene analysis."
        ) from exc
    return sc, sq


def parse_genes(value: str) -> list[str] | None:
    genes = [x.strip() for x in value.split(",") if x.strip()]
    return genes or None


def expression_vector(adata, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


def plot_gene(xy: np.ndarray, values: np.ndarray, gene: str, out: Path, point_size: float, invert_y: bool) -> None:
    order = np.argsort(values)
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    sca = ax.scatter(xy[order, 0], xy[order, 1], c=values[order], s=point_size, cmap="magma", linewidths=0)
    ax.set_aspect("equal")
    ax.axis("off")
    if invert_y:
        ax.invert_yaxis()
    ax.set_title(gene, fontsize=11, pad=6)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--genes", default="", help="Optional comma-separated genes. Default: highly variable genes.")
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=6)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sc, sq = require_packages()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    if args.spatial_key != "spatial":
        adata.obsm["spatial"] = adata.obsm[args.spatial_key].copy()

    genes = parse_genes(args.genes)
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    if genes is None:
        sc.pp.highly_variable_genes(work, n_top_genes=min(args.n_hvg, work.n_vars), flavor="seurat")
        genes = work.var_names[work.var.get("highly_variable", pd.Series(False, index=work.var_names)).to_numpy()].tolist()
    genes = [gene for gene in genes if gene in work.var_names]
    if not genes:
        raise SystemExit("No requested genes were found in the AnnData object.")

    sq.gr.spatial_neighbors(work, coord_type="generic", n_neighs=args.n_neighbors)
    sq.gr.spatial_autocorr(work, mode="moran", genes=genes)
    moran = work.uns.get("moranI")
    if moran is None or len(moran) == 0:
        raise SystemExit("squidpy did not return `adata.uns['moranI']`.")
    moran.to_csv(outdir / "spatially_variable_genes_moran.tsv", sep="\t")
    top = moran.sort_values("I", ascending=False).head(args.top_n)
    top.to_csv(outdir / "top_spatially_variable_genes.tsv", sep="\t")

    xy = np.asarray(work.obsm["spatial"])
    for gene in top.index:
        safe = gene.replace("/", "_").replace(" ", "_")
        plot_gene(xy, expression_vector(work, gene), gene, outdir / f"{safe}_spatial_expression.pdf", args.point_size, args.invert_y)
    print(f"Wrote spatially variable gene outputs to {outdir}")


if __name__ == "__main__":
    main()
