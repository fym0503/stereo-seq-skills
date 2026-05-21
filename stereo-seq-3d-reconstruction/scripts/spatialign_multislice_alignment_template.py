#!/usr/bin/env python3
"""Run the published spatiAlign-style multi-slice alignment workflow.

Template provenance:
Paper: spatiAlign: an unsupervised contrastive learning model for data
integration and alignment of spatially resolved transcriptomics
Paper DOI: 10.1093/gigascience/giae042
Original code: STOmics/Spatialign, demo.py
Reusable success pattern: keep the alignment wrapper parameterized, save the
trained/aligned output directory, and avoid hiding training settings.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def require_spatialign():
    try:
        from spatialign import Spatialign
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package/module `spatialign`. Install the published "
            "Spatialign package/repository before running this template. "
            "Blocked step: multi-slice spatial alignment."
        ) from exc
    return Spatialign


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-paths", nargs="+", required=True, help="Input slice files accepted by Spatialign.")
    parser.add_argument("--save-path", required=True, help="Directory for Spatialign outputs.")
    parser.add_argument("--min-genes", type=int, default=20)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=100)
    parser.add_argument("--n-neigh", type=int, default=15)
    parser.add_argument("--latent-dims", type=int, default=100)
    parser.add_argument("--max-epoch", type=int, default=500)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--patient", type=int, default=15)
    parser.add_argument("--tau1", type=float, default=0.2)
    parser.add_argument("--tau2", type=float, default=1.0)
    parser.add_argument("--tau3", type=float, default=0.5)
    parser.add_argument("--gpu", default=0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    Spatialign = require_spatialign()
    save_path = Path(args.save_path)
    save_path.mkdir(parents=True, exist_ok=True)

    model = Spatialign(
        *args.data_paths,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        batch_key=args.batch_key,
        is_norm_log=True,
        is_scale=False,
        is_hvg=False,
        is_reduce=False,
        n_pcs=args.n_pcs,
        n_hvg=args.n_hvg,
        n_neigh=args.n_neigh,
        is_undirected=True,
        latent_dims=args.latent_dims,
        tau1=args.tau1,
        tau2=args.tau2,
        tau3=args.tau3,
        is_verbose=True,
        seed=args.seed,
        gpu=args.gpu,
        save_path=str(save_path),
    )
    model.train(lr=args.lr, max_epoch=args.max_epoch, alpha=args.alpha, patient=args.patient)
    model.alignment()
    print(f"Wrote Spatialign outputs to {save_path}")


if __name__ == "__main__":
    main()
