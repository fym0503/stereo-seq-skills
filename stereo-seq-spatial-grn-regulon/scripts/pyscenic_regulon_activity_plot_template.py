#!/usr/bin/env python3
"""Plot pySCENIC/AUCell regulon activity in Stereo-seq coordinates.

Template provenance:
Stereo-seq corpus evidence: pySCENIC/SCENIC regulon analysis appears in 25
local Stereo-seq papers, including mouse placentation
(DOI: 10.1038/s41421-024-00740-6) and mouse organogenesis
(DOI: 10.1038/s41467-023-40155-7).
Representative code source: gpenglab/STcomm paper evidence plus pySCENIC tool
source, https://github.com/aertslab/pySCENIC
Reusable success pattern: reuse an existing AUCell matrix, select top
cluster-specific regulons, and export readable spatial maps and heatmaps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def require_scanpy():
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy`. Install it before running this "
            "template. Blocked step: loading AnnData for regulon activity plots."
        ) from exc
    return sc


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", index_col=0)


def top_regulons_by_group(auc: pd.DataFrame, groups: pd.Series, top_n: int) -> list[str]:
    grouped = auc.groupby(groups.astype(str)).mean()
    z = (grouped - grouped.mean(axis=0)) / grouped.std(axis=0).replace(0, np.nan)
    chosen: list[str] = []
    for _, row in z.iterrows():
        chosen.extend(row.sort_values(ascending=False).head(top_n).index.tolist())
    return list(dict.fromkeys(chosen))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="AnnData with obs labels and obsm spatial coordinates.")
    parser.add_argument("--auc", required=True, help="AUCell/regulon activity CSV/TSV, cells x regulons.")
    parser.add_argument("--cluster-key", default="annotation")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--outdir", default="pyscenic_regulon_plots")
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--point-size", type=float, default=2.0)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sns.set_theme(style="white", font="Arial")

    sc = require_scanpy()
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm spatial key: {args.spatial_key}")
    if args.cluster_key not in adata.obs:
        raise SystemExit(f"AnnData is missing obs cluster/domain key: {args.cluster_key}")

    auc = read_table(args.auc)
    common = adata.obs_names.intersection(auc.index)
    if len(common) == 0:
        raise SystemExit("No overlapping cells/bins between AnnData obs_names and AUCell matrix index.")
    adata = adata[common].copy()
    auc = auc.loc[common].copy()

    top_regs = top_regulons_by_group(auc, adata.obs[args.cluster_key], args.top_n)
    (outdir / "top_regulons.txt").write_text("\n".join(top_regs) + "\n", encoding="utf-8")

    grouped = auc[top_regs].groupby(adata.obs[args.cluster_key].astype(str)).mean()
    z = (grouped - grouped.mean(axis=0)) / grouped.std(axis=0).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(max(6, 0.36 * len(top_regs) + 2), max(4, 0.32 * len(grouped) + 2)))
    sns.heatmap(z.fillna(0), cmap="RdBu_r", center=0, linewidths=0.35, linecolor="white", ax=ax, cbar_kws={"label": "Regulon activity z-score"})
    ax.set_xlabel("Regulon")
    ax.set_ylabel(args.cluster_key)
    ax.tick_params(axis="x", labelrotation=45, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "regulon_activity_heatmap.pdf", dpi=300)
    plt.close(fig)

    xy = adata.obsm[args.spatial_key]
    for reg in top_regs:
        values = auc[reg]
        order = np.argsort(values.to_numpy())
        safe = reg.replace("/", "_").replace(" ", "_").replace("(+)", "")
        fig, ax = plt.subplots(figsize=(5.4, 5.2))
        sca = ax.scatter(xy[order, 0], xy[order, 1], c=values.to_numpy()[order], s=args.point_size, cmap="plasma", linewidths=0)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(reg, fontsize=11)
        fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
        fig.tight_layout()
        fig.savefig(outdir / f"{safe}_spatial_activity.pdf", dpi=300)
        plt.close(fig)

    print(f"Wrote pySCENIC/AUCell regulon plots to {outdir}")


if __name__ == "__main__":
    main()
