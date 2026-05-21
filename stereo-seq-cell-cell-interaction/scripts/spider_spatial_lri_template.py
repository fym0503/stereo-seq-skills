#!/usr/bin/env python3
"""Run a compact SPIDER spatial ligand-receptor interaction workflow.

Template provenance:
Paper: Finding spatially variable ligand-receptor interactions with functional
support from downstream genes
Paper DOI: 10.1038/s41467-025-62988-0
Original code: deepomicslab/SPIDER README.md and SPIDER-paper notebooks
Reusable success pattern: construct spatial interaction interfaces, score TF
support, test spatially variable LR interactions, and export pattern plots.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")


def require_spider():
    try:
        import anndata as ad
        import spider
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `spider-st`/`spider`. Install the SPIDER "
            "environment before running this template. Blocked step: spatially "
            "variable ligand-receptor interaction analysis."
        ) from exc
    return ad, spider


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Input AnnData .h5ad with spatial coordinates and cluster labels.")
    parser.add_argument("--cluster-key", default="celltype")
    parser.add_argument("--outdir", default="spider_lri_out")
    parser.add_argument("--species", choices=["human", "mouse"], default="mouse")
    parser.add_argument("--is-sc", action="store_true", help="Set if input is single-cell/cellbin scale.")
    parser.add_argument("--r-path", default="R", help="R executable path needed by SPIDER SVI tests.")
    parser.add_argument("--threshold", type=float, default=0.3)
    parser.add_argument("--alpha", type=float, default=0.3)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--show-svi", type=int, default=10)
    parser.add_argument("--spot-size", type=float, default=10)
    args = parser.parse_args()

    ad, spider = require_spider()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if args.cluster_key not in adata.obs:
        raise SystemExit(f"AnnData is missing obs cluster key: {args.cluster_key}")
    if "spatial" not in adata.obsm:
        raise SystemExit("AnnData is missing obsm['spatial'].")

    adata.uns["cluster_key"] = args.cluster_key
    adata.uns["is_human"] = args.species == "human"
    adata.uns["is_sc"] = args.is_sc

    op = spider.SPIDER()
    idata = op.prep(
        adata,
        cluster_key=args.cluster_key,
        is_human=adata.uns["is_human"],
        is_sc=adata.uns["is_sc"],
        itermax=1000,
        imputation=True,
        normalize_total=True,
    )
    op.svi.tf_corr(idata, adata, adata.uns["is_human"], str(outdir), threshold=args.threshold)
    idata, meta_idata = op.find_svi(idata, str(outdir), args.r_path, alpha=args.alpha, overwrite=True, n_jobs=args.n_jobs, svi_number=0)
    op.vis.pattern_LRI(idata, show_SVI=args.show_svi, spot_size=args.spot_size)
    svi_df, svi_df_strict = op.svi.combine_SVI(idata, threshold=0.01)
    svi_df.to_csv(outdir / "spider_svi_all.csv")
    svi_df_strict.to_csv(outdir / "spider_svi_strict.csv")
    idata.write_h5ad(outdir / "spider_interaction_data.h5ad")
    meta_idata.write_h5ad(outdir / "spider_meta_interaction_data.h5ad")
    print(f"Wrote SPIDER outputs to {outdir}")


if __name__ == "__main__":
    main()
