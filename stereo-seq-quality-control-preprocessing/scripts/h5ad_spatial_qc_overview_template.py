#!/usr/bin/env python3
"""Stereo-seq h5ad spatial QC overview template.

Template provenance:
Human endometrium PCOS atlas, DOI 10.1038/s41591-025-03592-z,
https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R, original file
`Spatial RNA-seq analysis/00.1_Stereopy_binning.py`.
Human cortex single-cell resolution atlas, DOI 10.1038/s41467-025-62793-9,
https://github.com/lcy1364/Cortex-Atlas-Code, original files
`src/STEREO/GEM/stat.R` and `src/STEREO/1_cellbin/1_gemToH5ad.py`.
GF/SPF cecum atlas, DOI 10.1016/j.isci.2024.108941,
https://github.com/1014723815/GF_SPF_cecum, original file
`Spatial transcriptome.R`.

Reusable success pattern: calculate count/feature QC before interpretation,
plot the metrics in tissue coordinates with equal aspect, keep legends/colorbars
outside the tissue body, and export manuscript-friendly PDF outputs.
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
            "Missing Python package `scanpy`. Install or activate an environment "
            "with scanpy before reading .h5ad input. Blocked step: h5ad spatial QC."
        ) from exc
    return sc


def matrix_sum(x, axis: int) -> np.ndarray:
    values = x.sum(axis=axis)
    return np.asarray(values).ravel()


def detected_genes(x) -> np.ndarray:
    values = (x > 0).sum(axis=1)
    return np.asarray(values).ravel()


def get_gene_vector(adata, gene_mask: np.ndarray) -> np.ndarray | None:
    if gene_mask.sum() == 0:
        return None
    values = adata.X[:, gene_mask].sum(axis=1)
    return np.asarray(values).ravel()


def split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def plot_hist(values: pd.Series, out: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.hist(values.dropna(), bins=60, color="#4E79A7", edgecolor="white", linewidth=0.25)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Spatial units", fontsize=10)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def plot_continuous_spatial(
    xy: np.ndarray,
    values: pd.Series,
    out: Path,
    title: str,
    point_size: float,
    invert_y: bool,
) -> None:
    vec = pd.to_numeric(values, errors="coerce").to_numpy()
    finite = np.isfinite(vec)
    if finite.sum() == 0:
        return
    lo, hi = np.nanquantile(vec[finite], [0.01, 0.99])
    clipped = np.clip(vec, lo, hi)
    order = np.argsort(clipped)
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    sca = ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=clipped[order],
        s=point_size,
        cmap="viridis",
        linewidths=0,
        alpha=0.95,
    )
    ax.set_aspect("equal")
    ax.axis("off")
    if invert_y:
        ax.invert_yaxis()
    ax.set_title(title, fontsize=11, pad=6)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def plot_categorical_spatial(
    xy: np.ndarray,
    labels: pd.Series,
    out: Path,
    title: str,
    point_size: float,
    invert_y: bool,
    max_categories: int,
) -> None:
    values = labels.astype(str).fillna("NA")
    counts = values.value_counts()
    keep = counts.head(max_categories).index
    values = values.where(values.isin(keep), other="other")
    categories = list(values.value_counts().index)
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(categories), 1)))
    palette = dict(zip(categories, colors))
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for category in categories:
        idx = values == category
        ax.scatter(
            xy[idx.to_numpy(), 0],
            xy[idx.to_numpy(), 1],
            s=point_size,
            color=palette[category],
            label=category,
            linewidths=0,
            alpha=0.95,
        )
    ax.set_aspect("equal")
    ax.axis("off")
    if invert_y:
        ax.invert_yaxis()
    ax.set_title(title, fontsize=11, pad=6)
    ax.legend(
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--obs-cols", default="", help="Comma-separated obs columns to map.")
    parser.add_argument("--prefix", default="spatial_qc")
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--max-categories", type=int, default=30)
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sc = require_scanpy()
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    xy = np.asarray(adata.obsm[args.spatial_key])
    if xy.shape[1] < 2:
        raise SystemExit(f"obsm['{args.spatial_key}'] must have at least two columns.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs.copy()
    obs["total_counts_template"] = matrix_sum(adata.X, axis=1)
    obs["detected_genes_template"] = detected_genes(adata.X)
    mt_mask = np.array([name.upper().startswith("MT-") for name in adata.var_names])
    mt_counts = get_gene_vector(adata, mt_mask)
    if mt_counts is not None:
        obs["mt_fraction_template"] = mt_counts / np.maximum(obs["total_counts_template"].to_numpy(), 1)

    obs[["total_counts_template", "detected_genes_template"] + (["mt_fraction_template"] if "mt_fraction_template" in obs else [])].to_csv(
        outdir / f"{args.prefix}_obs_qc.tsv",
        sep="\t",
    )
    summary = {
        "n_spatial_units": adata.n_obs,
        "n_genes": adata.n_vars,
        "x_min": float(np.nanmin(xy[:, 0])),
        "x_max": float(np.nanmax(xy[:, 0])),
        "y_min": float(np.nanmin(xy[:, 1])),
        "y_max": float(np.nanmax(xy[:, 1])),
        "template_sources": (
            "Endo.R DOI 10.1038/s41591-025-03592-z; "
            "Cortex-Atlas-Code DOI 10.1038/s41467-025-62793-9; "
            "GF_SPF_cecum DOI 10.1016/j.isci.2024.108941"
        ),
    }
    pd.Series(summary).to_csv(outdir / f"{args.prefix}_summary.tsv", sep="\t", header=False)

    for metric in ["total_counts_template", "detected_genes_template", "mt_fraction_template"]:
        if metric not in obs:
            continue
        plot_hist(obs[metric], outdir / f"{args.prefix}_{metric}_hist.pdf", metric)
        plot_continuous_spatial(
            xy,
            obs[metric],
            outdir / f"{args.prefix}_{metric}_spatial.pdf",
            metric,
            args.point_size,
            args.invert_y,
        )

    for col in split_csv(args.obs_cols):
        if col not in obs:
            print(f"Skip missing obs column: {col}")
            continue
        if pd.api.types.is_numeric_dtype(obs[col]):
            plot_continuous_spatial(xy, obs[col], outdir / f"{args.prefix}_{col}_spatial.pdf", col, args.point_size, args.invert_y)
        else:
            plot_categorical_spatial(
                xy,
                obs[col],
                outdir / f"{args.prefix}_{col}_spatial.pdf",
                col,
                args.point_size,
                args.invert_y,
                args.max_categories,
            )
    print(f"Wrote h5ad spatial QC overview to {outdir}")


if __name__ == "__main__":
    main()
