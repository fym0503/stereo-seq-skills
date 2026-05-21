#!/usr/bin/env python3
"""Scanpy spatial domain and marker template for Stereo-seq h5ad objects.

Template provenance:
P09 AD prefrontal cortex, DOI 10.1038/s41467-024-54715-y,
Zenodo 10.5281/zenodo.14048103, original files `1_bin110_analysis.ipynb`
and `2_layer_annotation.ipynb`.
Human endometrium PCOS atlas, DOI 10.1038/s41591-025-03592-z,
https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R, original file
`Spatial RNA-seq analysis/00.1_Stereopy_binning.py`.
Human cortex single-cell resolution atlas, DOI 10.1038/s41467-025-62793-9,
https://github.com/lcy1364/Cortex-Atlas-Code, original domain files under
`src/STEREO/3_domainProcess/`.

Reusable success pattern: run an auditable unsupervised domain pass, export
domain labels and top markers together, and validate domains with equal-aspect
spatial maps plus a readable marker heatmap.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def require_scanpy():
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy`. Install or activate scanpy before "
            "running this domain template. Blocked step: unsupervised spatial domain analysis."
        ) from exc
    return sc


def stable_palette(labels: list[str]) -> dict[str, tuple[float, float, float, float]]:
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(labels), 1)))
    return dict(zip(labels, colors))


def plot_domains(xy: np.ndarray, labels: pd.Series, out: Path, point_size: float, invert_y: bool) -> None:
    categories = sorted(labels.astype(str).unique())
    palette = stable_palette(categories)
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for category in categories:
        idx = labels.astype(str) == category
        ax.scatter(xy[idx.to_numpy(), 0], xy[idx.to_numpy(), 1], s=point_size, color=palette[category], label=category, linewidths=0)
    ax.set_aspect("equal")
    ax.axis("off")
    if invert_y:
        ax.invert_yaxis()
    ax.legend(
        title="Domain",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
        title_fontsize=10,
        markerscale=max(2.5, 4 / max(point_size, 0.5)),
    )
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def marker_table(adata, group_key: str, top_n: int) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    records: list[dict[str, object]] = []
    for group in groups:
        names = result["names"][group][:top_n]
        scores = result["scores"][group][:top_n] if "scores" in result else [np.nan] * len(names)
        logfc = result["logfoldchanges"][group][:top_n] if "logfoldchanges" in result else [np.nan] * len(names)
        pvals = result["pvals_adj"][group][:top_n] if "pvals_adj" in result else [np.nan] * len(names)
        for rank, (gene, score, fc, pval) in enumerate(zip(names, scores, logfc, pvals), start=1):
            records.append(
                {
                    "group_key": group_key,
                    "domain": group,
                    "rank": rank,
                    "gene": gene,
                    "score": score,
                    "logfoldchange": fc,
                    "p_adj": pval,
                }
            )
    return pd.DataFrame(records)


def plot_marker_heatmap(adata, markers: list[str], group_key: str, out: Path) -> None:
    markers = [gene for gene in markers if gene in adata.var_names]
    if not markers:
        return
    x = adata[:, markers].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    expr = pd.DataFrame(np.asarray(x), index=adata.obs[group_key].astype(str), columns=markers)
    grouped = expr.groupby(level=0).mean()
    z = (grouped - grouped.mean(axis=0)) / grouped.std(axis=0).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(max(6, 0.28 * len(markers) + 2), max(4, 0.32 * len(grouped) + 1.5)))
    im = ax.imshow(z.fillna(0).to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(markers)))
    ax.set_xticklabels(markers, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(grouped.index)))
    ax.set_yticklabels(grouped.index, fontsize=9)
    ax.set_xlabel("Top marker genes", fontsize=10)
    ax.set_ylabel(group_key, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Mean expression z-score", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--cluster-key", default="template_domain")
    parser.add_argument("--use-existing-key", default="", help="Use this obs key instead of running Leiden.")
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--top-markers", type=int, default=8)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sc = require_scanpy()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")

    if args.use_existing_key:
        if args.use_existing_key not in adata.obs:
            raise SystemExit(f"AnnData is missing obs['{args.use_existing_key}'].")
        args.cluster_key = args.use_existing_key
    else:
        try:
            work = adata.copy()
            sc.pp.normalize_total(work, target_sum=1e4)
            sc.pp.log1p(work)
            sc.pp.highly_variable_genes(work, n_top_genes=min(3000, work.n_vars), flavor="seurat")
            if "highly_variable" in work.var:
                work = work[:, work.var["highly_variable"]].copy()
            sc.pp.scale(work, max_value=10)
            sc.tl.pca(work, n_comps=min(args.n_pcs, work.n_vars - 1, work.n_obs - 1))
            sc.pp.neighbors(work, n_neighbors=args.n_neighbors, n_pcs=min(args.n_pcs, work.obsm["X_pca"].shape[1]))
            sc.tl.leiden(work, resolution=args.resolution, key_added=args.cluster_key)
            adata.obs[args.cluster_key] = work.obs[args.cluster_key].astype(str).to_numpy()
        except Exception as exc:
            raise SystemExit(
                "Leiden/Scanpy domain discovery failed. Check that `leidenalg`, "
                "`igraph`, and scanpy dependencies are installed. Blocked step: "
                f"unsupervised domain clustering. Original error: {exc}"
            ) from exc

    sc.tl.rank_genes_groups(adata, groupby=args.cluster_key, method="wilcoxon")
    markers = marker_table(adata, args.cluster_key, args.top_markers)
    markers.to_csv(outdir / "domain_top_markers.tsv", sep="\t", index=False)
    adata.obs[[args.cluster_key]].to_csv(outdir / "domain_labels.tsv", sep="\t")
    adata.write_h5ad(outdir / "domain_result.h5ad")

    xy = np.asarray(adata.obsm[args.spatial_key])
    plot_domains(xy, adata.obs[args.cluster_key].astype(str), outdir / "domain_spatial_map.pdf", args.point_size, args.invert_y)
    heatmap_genes = list(dict.fromkeys(markers["gene"].head(args.top_markers * adata.obs[args.cluster_key].nunique()).astype(str)))
    plot_marker_heatmap(adata, heatmap_genes, args.cluster_key, outdir / "domain_marker_heatmap.pdf")
    print(f"Wrote scanpy spatial domain outputs to {outdir}")


if __name__ == "__main__":
    main()
