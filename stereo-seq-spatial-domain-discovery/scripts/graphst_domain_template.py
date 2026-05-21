#!/usr/bin/env python3
"""GraphST spatial-domain discovery template.

Template provenance:
Paper: Spatially informed clustering, integration, and deconvolution of spatial
transcriptomics with GraphST
Paper DOI: 10.1038/s41467-023-36796-3
Original code: https://github.com/JinmiaoChenLab/GraphST,
`GraphST/GraphST.py` and `GraphST/utils.py`.

Reusable success pattern: learn a spatial graph representation from gene
expression and coordinates, cluster the learned embedding into spatial domains,
optionally refine labels by spatial neighbors, and export equal-aspect domain
maps.
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


def require_graphst():
    try:
        import torch
        from GraphST import GraphST
        from GraphST.utils import clustering
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `GraphST` and/or PyTorch. Install or activate "
            "the GraphST environment before running. Blocked step: GraphST domain discovery."
        ) from exc
    return torch, GraphST, clustering


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", default="graphst_domain_out")
    parser.add_argument("--n-clusters", type=int, required=True)
    parser.add_argument("--method", choices=["mclust", "leiden", "louvain"], default="mclust")
    parser.add_argument("--datatype", default="Stereo", choices=["Stereo", "Slide", "10X"])
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--refinement", action="store_true")
    parser.add_argument("--radius", type=int, default=50)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--point-size", type=float, default=1.0)
    args = parser.parse_args()

    torch, GraphST, clustering = require_graphst()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    if args.spatial_key != "spatial":
        adata.obsm["spatial"] = adata.obsm[args.spatial_key].copy()

    if args.device == "auto":
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    model = GraphST.GraphST(adata, device=device, epochs=args.epochs, datatype=args.datatype)
    adata = model.train()
    clustering(
        adata,
        n_clusters=args.n_clusters,
        radius=args.radius,
        key="emb",
        method=args.method,
        refinement=args.refinement,
    )
    adata.write_h5ad(outdir / "graphst_domain_result.h5ad")
    adata.obs[["domain"]].to_csv(outdir / "graphst_domain_labels.csv")

    xy = adata.obsm["spatial"]
    labels = adata.obs["domain"].astype(str)
    palette = dict(zip(sorted(labels.unique()), plt.cm.tab20(np.linspace(0, 1, len(labels.unique())))))
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    for label in sorted(labels.unique()):
        idx = labels == label
        ax.scatter(xy[idx, 0], xy[idx, 1], s=args.point_size, color=palette[label], label=label, linewidths=0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(
        title="GraphST domain",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
        title_fontsize=10,
        markerscale=max(3, 4 / max(args.point_size, 0.5)),
    )
    fig.tight_layout()
    fig.savefig(outdir / "graphst_spatial_domains.pdf", dpi=300)
    plt.close(fig)
    print(f"Wrote GraphST domain outputs to {outdir}")


if __name__ == "__main__":
    main()
