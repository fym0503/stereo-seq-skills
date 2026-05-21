#!/usr/bin/env python3
"""Run a compact cell2location mapping workflow for Stereo-seq-like AnnData.

Template provenance:
Paper: Effects of flora deficiency on the structure and function of the large intestine
Paper DOI: 10.1016/j.isci.2024.108941
Original code: 1014723815/GF_SPF_cecum, cell2location.py
Reusable success pattern: train reference signatures, map cell abundance to
spatial locations, export posterior abundance matrix for downstream spatial plots.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def require_cell2location():
    try:
        import cell2location
        from cell2location.models import Cell2location, RegressionModel
        from cell2location.utils.filtering import filter_genes
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `cell2location`. Install it in the active "
            "environment before running this template. Blocked step: "
            "cell2location reference regression and spatial abundance mapping."
        ) from exc
    return cell2location, Cell2location, RegressionModel, filter_genes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial", required=True, help="Spatial AnnData .h5ad")
    parser.add_argument("--reference", required=True, help="sc/snRNA reference AnnData .h5ad")
    parser.add_argument("--label-key", default="celltype", help="Reference obs column with labels")
    parser.add_argument("--batch-key", default=None, help="Optional reference obs batch column")
    parser.add_argument("--outdir", default="cell2location_out")
    parser.add_argument("--expected-cells", type=float, default=4)
    parser.add_argument("--use-gpu", action="store_true")
    args = parser.parse_args()

    cell2location, Cell2location, RegressionModel, filter_genes = require_cell2location()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata_sp = sc.read_h5ad(args.spatial)
    adata_ref = sc.read_h5ad(args.reference)
    if args.label_key not in adata_ref.obs:
        raise SystemExit(f"Reference obs is missing label key: {args.label_key}")

    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_ref = adata_ref[:, selected].copy()

    RegressionModel.setup_anndata(adata=adata_ref, batch_key=args.batch_key, labels_key=args.label_key)
    ref_model = RegressionModel(adata_ref)
    ref_model.train(max_epochs=200, batch_size=2500, train_size=1, lr=0.005, use_gpu=args.use_gpu)
    adata_ref = ref_model.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500, "use_gpu": args.use_gpu},
    )
    ref_model.save(outdir / "reference_regression_model", overwrite=True)

    if "means_per_cluster_mu_fg" in adata_ref.varm:
        signature = adata_ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{name}" for name in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    else:
        signature = adata_ref.var[
            [f"means_per_cluster_mu_fg_{name}" for name in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    signature.columns = adata_ref.uns["mod"]["factor_names"]

    adata_sp.var["mt_gene"] = adata_sp.var_names.str.upper().str.startswith("MT-")
    adata_sp = adata_sp[:, ~adata_sp.var["mt_gene"].to_numpy()].copy()
    common = np.intersect1d(adata_sp.var_names, signature.index)
    adata_sp = adata_sp[:, common].copy()
    signature = signature.loc[common].copy()

    Cell2location.setup_anndata(adata=adata_sp)
    spatial_model = Cell2location(
        adata_sp,
        cell_state_df=signature,
        N_cells_per_location=args.expected_cells,
        detection_alpha=20,
    )
    spatial_model.train(max_epochs=8000, batch_size=None, train_size=1, use_gpu=args.use_gpu)
    adata_sp = spatial_model.export_posterior(
        adata_sp,
        sample_kwargs={"num_samples": 500, "batch_size": adata_sp.n_obs, "use_gpu": args.use_gpu},
    )
    spatial_model.save(outdir / "spatial_cell2location_model", overwrite=True)
    adata_sp.write_h5ad(outdir / "spatial_cell2location_mapped.h5ad")

    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key in adata_sp.obsm:
        abundance = adata_sp.obsm[abundance_key]
        if not hasattr(abundance, "to_csv"):
            factor_names = adata_sp.uns.get("mod", {}).get("factor_names", [f"factor_{i}" for i in range(abundance.shape[1])])
            abundance = pd.DataFrame(abundance, index=adata_sp.obs_names, columns=factor_names)
        abundance.to_csv(outdir / "cell_abundance_q05.csv")
        if "spatial" in adata_sp.obsm:
            top_labels = abundance.sum(axis=0).sort_values(ascending=False).head(6).index
            for label in top_labels:
                values = abundance[label].to_numpy()
                fig, ax = plt.subplots(figsize=(5.5, 5.5))
                sca = ax.scatter(
                    adata_sp.obsm["spatial"][:, 0],
                    adata_sp.obsm["spatial"][:, 1],
                    c=values,
                    s=2,
                    cmap="magma",
                    linewidths=0,
                )
                ax.set_aspect("equal")
                ax.axis("off")
                ax.set_title(str(label), fontsize=11, fontfamily="Arial")
                fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
                fig.tight_layout()
                fig.savefig(outdir / f"abundance_{label}.pdf", dpi=300)
                plt.close(fig)

    print(f"Wrote cell2location outputs to {outdir}")


if __name__ == "__main__":
    main()
