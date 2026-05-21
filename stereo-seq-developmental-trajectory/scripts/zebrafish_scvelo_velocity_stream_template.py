#!/usr/bin/env python3
"""Run a compact scVelo RNA velocity stream workflow.

Template provenance:
Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the
regenerating zebrafish heart
Paper DOI: 10.1038/s41467-025-59070-0
Original code: BGI-Qingdao/ZebrafishHeartRegeneration_project,
  06. RNA velocity analysis/03.velocity.py
Reusable success pattern: merge loom velocity layers with metadata/UMAP,
run scVelo, and export readable stream plots by cell state.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def require_scvelo():
    try:
        import scvelo as scv
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scvelo`. Install it before running this "
            "template. Blocked step: RNA velocity preprocessing and stream plots."
        ) from exc
    return scv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loom", required=True, help="Input loom file with spliced/unspliced layers.")
    parser.add_argument("--meta", required=True, help="CSV/TSV metadata table keyed by cell IDs.")
    parser.add_argument("--cell-id-col", default="CellID")
    parser.add_argument("--umap-cols", default="umap_1,umap_2", help="Two metadata columns for UMAP coordinates.")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--stage-col", default="", help="Optional stage/time column for an additional colored plot.")
    parser.add_argument("--outdir", default="velocity_out")
    parser.add_argument("--prefix", default="velocity")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--mode", choices=["stochastic", "dynamical"], default="dynamical")
    args = parser.parse_args()

    scv = require_scvelo()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams["font.family"] = "Arial"
    scv.settings.figdir = str(outdir)
    scv.settings.set_figure_params("scvelo", dpi=120, dpi_save=300, fontsize=10)

    adata = scv.read(args.loom)
    meta = pd.read_csv(args.meta, sep=None, engine="python")
    if args.cell_id_col not in meta:
        raise SystemExit(f"Metadata is missing cell ID column: {args.cell_id_col}")
    meta = meta.set_index(args.cell_id_col)

    common = adata.obs_names.intersection(meta.index)
    if len(common) == 0:
        raise SystemExit("No overlapping cell IDs between loom obs_names and metadata.")
    adata = adata[common].copy()
    adata.obs = adata.obs.join(meta.loc[common], how="left")

    umap_cols = [col.strip() for col in args.umap_cols.split(",")]
    if len(umap_cols) != 2 or any(col not in adata.obs for col in umap_cols):
        raise SystemExit(f"UMAP columns are missing from metadata: {umap_cols}")
    adata.obsm["X_umap"] = adata.obs[umap_cols].to_numpy()

    scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=args.n_top_genes)
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)

    adata.write(outdir / f"{args.prefix}_velocity.h5ad")
    color_keys = [args.celltype_col]
    if args.stage_col and args.stage_col in adata.obs:
        color_keys.append(args.stage_col)

    for color in color_keys:
        if color not in adata.obs:
            continue
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color=color,
            legend_loc="right margin",
            frameon=False,
            title=color,
            save=f"_{args.prefix}_{color}_stream.pdf",
            show=False,
        )
    scv.pl.scatter(
        adata,
        basis="umap",
        color="velocity_confidence",
        cmap="viridis",
        frameon=False,
        title="Velocity confidence",
        save=f"_{args.prefix}_confidence.pdf",
        show=False,
    )
    print(f"Wrote velocity outputs to {outdir}")


if __name__ == "__main__":
    main()
