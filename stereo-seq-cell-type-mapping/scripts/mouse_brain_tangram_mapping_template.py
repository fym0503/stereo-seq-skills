#!/usr/bin/env python3
"""Tangram cell/reference-to-space mapping template.

Template provenance:
Paper: Benchmarking mapping algorithms for cell-type annotating in mouse brain
by integrating single-nucleus RNA-seq and Stereo-seq data
Paper DOI: 10.1093/bib/bbae250
Original code: https://github.com/qyTao185/Benchmarking-Mapping-Algorithms,
`Code in python/tangram_circle.ipynb`.
Related Stereo-seq evidence: stTransfer paper, DOI 10.1016/j.crmeth.2025.101205,
lists Tangram among evaluated single-cell-to-spatial mapping methods.

Reusable success pattern: rank reference markers, intersect genes with spatial
data, run `tg.pp_adatas`, map cells to space with RNA-count density prior, then
export projected cell-type probabilities for downstream spatial plotting.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def require_tangram():
    try:
        import tangram as tg
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `tangram-sc`/`tangram`. Install or activate "
            "an environment containing Tangram before running. Blocked step: "
            "Tangram cell-to-space mapping."
        ) from exc
    return tg


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial", required=True, help="Spatial AnnData .h5ad")
    parser.add_argument("--reference", required=True, help="Reference sc/snRNA AnnData .h5ad")
    parser.add_argument("--label-key", default="annotation")
    parser.add_argument("--outdir", default="tangram_mapping_out")
    parser.add_argument("--mode", default="cells", choices=["cells", "clusters"])
    parser.add_argument("--top-marker-n", type=int, default=100)
    parser.add_argument("--density-prior", default="rna_count_based")
    parser.add_argument("--num-epochs", type=int, default=100)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--point-size", type=float, default=1.2)
    args = parser.parse_args()

    tg = require_tangram()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

    adata_sc = sc.read_h5ad(args.reference)
    adata_st = sc.read_h5ad(args.spatial)
    if args.label_key not in adata_sc.obs:
        raise SystemExit(f"Reference obs is missing label key: {args.label_key}")

    start = time.time()
    sc.tl.rank_genes_groups(adata_sc, groupby=args.label_key, use_raw=False)
    marker_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[: args.top_marker_n, :]
    markers = list(pd.unique(marker_df.melt().value.dropna()))
    common_markers = sorted(set(markers).intersection(adata_sc.var_names).intersection(adata_st.var_names))
    if not common_markers:
        raise SystemExit("No marker genes overlap between reference and spatial AnnData.")

    tg.pp_adatas(adata_sc, adata_st, genes=common_markers)
    map_kwargs = {
        "mode": args.mode,
        "density_prior": args.density_prior,
        "num_epochs": args.num_epochs,
        "device": args.device,
    }
    if args.mode == "clusters":
        map_kwargs["cluster_label"] = args.label_key
    ad_map = tg.map_cells_to_space(adata_sc, adata_st, **map_kwargs)
    tg.project_cell_annotations(ad_map, adata_st, annotation=args.label_key)

    pred = adata_st.obsm["tangram_ct_pred"].copy()
    if not isinstance(pred, pd.DataFrame):
        pred = pd.DataFrame(pred, index=adata_st.obs_names)
    pred.to_csv(outdir / "tangram_celltype_prediction.csv")
    pd.Series({"seconds": time.time() - start, "n_markers": len(common_markers)}).to_csv(outdir / "tangram_runtime.csv")
    adata_st.write_h5ad(outdir / "spatial_tangram_projected.h5ad")

    if args.spatial_key in adata_st.obsm:
        xy = adata_st.obsm[args.spatial_key]
        top_labels = pred.sum(axis=0).sort_values(ascending=False).head(6).index
        for label in top_labels:
            values = pred[label].to_numpy()
            fig, ax = plt.subplots(figsize=(5.4, 5.2))
            sca = ax.scatter(xy[:, 0], xy[:, 1], c=values, s=args.point_size, cmap="magma", linewidths=0)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(str(label), fontsize=11)
            fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
            fig.tight_layout()
            safe = str(label).replace("/", "_").replace(" ", "_")
            fig.savefig(outdir / f"tangram_{safe}_spatial_probability.pdf", dpi=300)
            plt.close(fig)

    print(f"Wrote Tangram mapping outputs to {outdir}")


if __name__ == "__main__":
    main()
