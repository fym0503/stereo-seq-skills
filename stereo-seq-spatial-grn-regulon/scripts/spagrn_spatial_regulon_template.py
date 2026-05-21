#!/usr/bin/env python3
"""Run a compact SpaGRN workflow and export regulon plots.

Template provenance:
Paper: SpaGRN: Investigating spatially informed regulatory paths for spatially
resolved transcriptomics data
Paper DOI: 10.1016/j.cels.2025.101243
Original code: BGI-Qingdao/SpaGRN,
  docs/source/content/01_Basic_Usage.rst and
  docs/source/Tutorials/stereo_seq_mouse_brain_*.ipynb
Reusable success pattern: infer spatially informed TF regulons from AnnData
using spatial coordinates, then plot top regulon activity in tissue space.
"""

from __future__ import annotations

import argparse
import json
from multiprocessing import cpu_count
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def require_spagrn():
    try:
        import scanpy as sc
        import spagrn.plot as prn
        from spagrn.regulatory_network import InferNetwork as irn
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python packages `spagrn` and/or `scanpy`. Install SpaGRN "
            "and its pySCENIC dependencies before running this template. "
            "Blocked step: spatial GRN inference and regulon plotting."
        ) from exc
    return sc, prn, irn


def auc_dataframe(adata) -> pd.DataFrame:
    auc_raw = adata.obsm["auc_mtx"]
    if isinstance(auc_raw, pd.DataFrame):
        return auc_raw.copy()
    columns = list(adata.uns.get("regulon_dict", {}).keys())
    if len(columns) != auc_raw.shape[1]:
        columns = [f"regulon_{i}" for i in range(auc_raw.shape[1])]
    return pd.DataFrame(auc_raw, index=adata.obs_names, columns=columns)


def choose_top_regulons(adata, cluster_key: str, top_n: int) -> list[str]:
    auc = auc_dataframe(adata)
    if cluster_key not in adata.obs:
        return auc.mean(axis=0).sort_values(ascending=False).head(top_n).index.tolist()
    grouped = auc.groupby(adata.obs[cluster_key].astype(str)).mean()
    z = (grouped - grouped.mean(axis=0)) / grouped.std(axis=0).replace(0, np.nan)
    ranked = []
    for _, row in z.iterrows():
        ranked.extend(row.sort_values(ascending=False).head(top_n).index.tolist())
    return list(dict.fromkeys(ranked))[: max(top_n, 1) * min(len(grouped), 4)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Input AnnData .h5ad with spatial coordinates.")
    parser.add_argument("--database", required=True, help="cisTarget ranking database .feather.")
    parser.add_argument("--motif", required=True, help="Motif annotation .tbl.")
    parser.add_argument("--tfs", required=True, help="TF list text file.")
    parser.add_argument("--cluster-key", default="annotation")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--layer-key", default="raw_counts")
    parser.add_argument("--project-name", default="stereo_spagrn")
    parser.add_argument("--outdir", default="spagrn_out")
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--model", choices=["bernoulli", "danb", "normal", "none"], default="bernoulli")
    parser.add_argument("--mode", choices=["moran", "geary"], default="moran")
    parser.add_argument("--methods", default="FDR_I,FDR_C,FDR_G")
    parser.add_argument("--operation", choices=["intersection", "union"], default="intersection")
    parser.add_argument("--num-workers", type=int, default=min(cpu_count(), 8))
    parser.add_argument("--rank-threshold", type=int, default=9000)
    parser.add_argument("--prune-auc-threshold", type=float, default=0.05)
    parser.add_argument("--auc-threshold", type=float, default=0.05)
    parser.add_argument("--top-regulons", type=int, default=6)
    args = parser.parse_args()

    _, prn, irn = require_spagrn()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

    data = irn.read_file(args.h5ad)
    if args.spatial_key not in data.obsm:
        raise SystemExit(f"AnnData is missing obsm spatial key: {args.spatial_key}")
    if args.cluster_key not in data.obs:
        raise SystemExit(f"AnnData is missing obs cluster/domain key: {args.cluster_key}")

    data = irn.preprocess(data)
    grn = irn(data, project_name=args.project_name)
    grn.add_params(
        {
            "prune_auc_threshold": args.prune_auc_threshold,
            "rank_threshold": args.rank_threshold,
            "auc_threshold": args.auc_threshold,
        }
    )
    methods = [method.strip() for method in args.methods.split(",") if method.strip()]
    grn.infer(
        args.database,
        args.motif,
        args.tfs,
        gene_list=None,
        num_workers=args.num_workers,
        cache=False,
        output_dir=str(outdir),
        save_tmp=True,
        layers=args.layer_key,
        latent_obsm_key=args.spatial_key,
        model=args.model,
        n_neighbors=args.n_neighbors,
        methods=methods,
        operation=args.operation,
        mode=args.mode,
        cluster_label=args.cluster_key,
    )

    result_h5ad = outdir / f"{args.project_name}_spagrn.h5ad"
    adata = irn.read_file(str(result_h5ad))
    auc_mtx = auc_dataframe(adata)
    auc_mtx.to_csv(outdir / f"{args.project_name}_auc_mtx.csv")
    (outdir / f"{args.project_name}_regulons.json").write_text(
        json.dumps(adata.uns.get("regulon_dict", {}), indent=2, default=list),
        encoding="utf-8",
    )

    top_regulons = choose_top_regulons(adata, args.cluster_key, args.top_regulons)
    for regulon in top_regulons:
        safe = regulon.replace("/", "_").replace(" ", "_").replace("(+)", "")
        prn.plot_2d_reg(
            adata,
            auc_mtx,
            reg_name=regulon,
            pos_label=args.spatial_key,
            fn=str(outdir / f"{args.project_name}_{safe}_spatial.pdf"),
            s=2,
        )

    try:
        prn.auc_heatmap(
            adata,
            auc_mtx,
            cluster_label=args.cluster_key,
            topn=args.top_regulons,
            subset=False,
            save=True,
            fn=str(outdir / f"{args.project_name}_regulon_heatmap.pdf"),
            legend_fn=str(outdir / f"{args.project_name}_regulon_legend.pdf"),
            figsize=(8, 7),
        )
    except Exception as exc:
        print(f"Warning: SpaGRN heatmap plotting failed after inference: {exc}")

    print(f"Wrote SpaGRN outputs to {outdir}")


if __name__ == "__main__":
    main()
