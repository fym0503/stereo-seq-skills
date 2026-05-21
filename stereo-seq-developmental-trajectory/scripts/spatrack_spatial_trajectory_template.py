#!/usr/bin/env python3
"""spaTrack spatial trajectory template.

Template provenance:
Paper: Inferring cell trajectories of spatial transcriptomics via optimal
transport analysis
Paper DOI: 10.1016/j.cels.2025.101194
Code DOI: 10.5281/zenodo.14214597
Original code: https://github.com/yzf072/spaTrack,
`docs/source/notebooks/04.ST_data_of_mouse midbrain_with_multiple_times.ipynb`
and `spaTrack/single_time/velocity.py`.

Reusable success pattern: estimate an optimal-transport transition matrix from
gene expression and spatial distance, choose explicit starting cells, calculate
pseudotime, and plot equal-aspect spatial pseudotime/trajectory stream panels.
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


def require_spatrack():
    try:
        import spaTrack as spt
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `spaTrack`. Install or activate the spaTrack "
            "environment before running. Blocked step: spatial OT trajectory."
        ) from exc
    return spt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", default="spatrack_trajectory_out")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--cluster-key", default="cluster")
    parser.add_argument("--start-cell-type", default="")
    parser.add_argument("--start-x", type=float, default=None)
    parser.add_argument("--start-y", type=float, default=None)
    parser.add_argument("--alpha-gene", type=float, default=0.5)
    parser.add_argument("--alpha-spatial", type=float, default=0.5)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neigh-pos", type=int, default=10)
    parser.add_argument("--grid-num", type=int, default=50)
    parser.add_argument("--point-size", type=float, default=1.2)
    args = parser.parse_args()

    spt = require_spatrack()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    adata.obsm["X_spatial"] = adata.obsm[args.spatial_key].copy()
    if args.cluster_key in adata.obs and args.cluster_key != "cluster":
        adata.obs["cluster"] = adata.obs[args.cluster_key].astype(str)
    elif "cluster" not in adata.obs:
        raise SystemExit("Provide --cluster-key or include obs['cluster'] for start-cell selection.")

    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.pca(adata, n_comps=max(args.n_pcs, 2), svd_solver="arpack")

    trans = spt.get_ot_matrix(
        adata,
        data_type="spatial",
        alpha1=args.alpha_gene,
        alpha2=args.alpha_spatial,
        n_pcs=args.n_pcs,
    )
    adata.obsp["trans"] = trans

    if args.start_cell_type:
        start_cells = spt.set_start_cells(adata, select_way="cell_type", cell_type=args.start_cell_type, basis="spatial")
    else:
        if args.start_x is None or args.start_y is None:
            raise SystemExit("Provide either --start-cell-type or both --start-x/--start-y.")
        start_cells = spt.set_start_cells(
            adata,
            select_way="coordinates",
            start_point=(args.start_x, args.start_y),
            basis="spatial",
        )
    if not start_cells:
        raise SystemExit("spaTrack start-cell selection returned no cells; adjust start coordinates or cell type.")

    adata.obs["ptime"] = spt.get_ptime(adata, start_cells=start_cells)
    p_grid, v_grid = spt.get_velocity(
        adata,
        basis="spatial",
        n_neigh_pos=args.n_neigh_pos,
        n_neigh_gene=0,
        grid_num=args.grid_num,
    )
    adata.write_h5ad(outdir / "spatrack_trajectory_result.h5ad")
    adata.obs[["cluster", "ptime"]].to_csv(outdir / "spatrack_pseudotime.csv")
    pd.DataFrame(trans, index=adata.obs_names, columns=adata.obs_names).to_csv(outdir / "spatrack_transition_matrix.csv")

    xy = adata.obsm["X_spatial"]
    order = np.argsort(adata.obs["ptime"].to_numpy())
    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    sca = ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=adata.obs["ptime"].to_numpy()[order],
        s=args.point_size,
        cmap="viridis",
        linewidths=0,
    )
    ax.scatter(xy[start_cells, 0], xy[start_cells, 1], s=18, c="black", marker="x", linewidths=0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("spaTrack pseudotime", fontsize=11)
    fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
    fig.tight_layout()
    fig.savefig(outdir / "spatrack_spatial_pseudotime.pdf", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    ax.scatter(xy[:, 0], xy[:, 1], c=adata.obs["ptime"], s=args.point_size, cmap="Greys", alpha=0.45, linewidths=0)
    ax.streamplot(
        p_grid[0],
        p_grid[1],
        v_grid[0],
        v_grid[1],
        color="#D62728",
        density=1.0,
        linewidth=0.7,
        arrowsize=0.8,
    )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("spaTrack spatial trajectory", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "spatrack_velocity_stream.pdf", dpi=300)
    plt.close(fig)

    print(f"Wrote spaTrack trajectory outputs to {outdir}")


if __name__ == "__main__":
    main()
