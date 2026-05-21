#!/usr/bin/env python3
"""Spatial pseudotime/state-transition plotting template.

Template provenance:
spaTrack, DOI 10.1016/j.cels.2025.101194, code DOI 10.5281/zenodo.14214597,
https://github.com/yzf072/spaTrack, original notebooks under
`docs/source/notebooks/04.ST_data_of_mouse midbrain_with_multiple_times.ipynb`.
ONTraC, DOI 10.1186/s13059-025-03588-5, code DOI 10.5281/zenodo.14171604,
https://github.com/gyuanlab/ONTraC and ONTraC_paper Stereo-seq midbrain
post-analysis notebooks.
Ascidian endostyle, DOI 10.1126/sciadv.adi9035,
https://github.com/lskfs/ascidian-endostyle, original files
`11.trajectory_projection/trajectory_mapping.py` and `monocle2.R`.

Reusable success pattern: plot externally validated pseudotime or state scores
in tissue coordinates, avoid claiming direction without time/velocity/start-cell
evidence, and export both a spatial map and per-state summaries.
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
            "Missing Python package `scanpy`. Install it before reading h5ad "
            "for spatial pseudotime plotting. Blocked figure: pseudotime map."
        ) from exc
    return sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--pseudotime-key", required=True)
    parser.add_argument("--state-key", default="")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--invert-y", action="store_true")
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    sc = require_scanpy()
    adata = sc.read_h5ad(args.h5ad)
    if args.spatial_key not in adata.obsm:
        raise SystemExit(f"AnnData is missing obsm['{args.spatial_key}'].")
    if args.pseudotime_key not in adata.obs:
        raise SystemExit(f"AnnData is missing obs['{args.pseudotime_key}'].")
    xy = np.asarray(adata.obsm[args.spatial_key])
    pt = pd.to_numeric(adata.obs[args.pseudotime_key], errors="coerce")
    order = np.argsort(pt.fillna(pt.min()).to_numpy())
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    sca = ax.scatter(xy[order, 0], xy[order, 1], c=pt.to_numpy()[order], s=args.point_size, cmap="viridis", linewidths=0)
    ax.set_aspect("equal")
    ax.axis("off")
    if args.invert_y:
        ax.invert_yaxis()
    ax.set_title(args.pseudotime_key, fontsize=11, pad=6)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Pseudotime", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_spatial_pseudotime.pdf", dpi=300)
    plt.close(fig)

    table = pd.DataFrame({"pseudotime": pt.to_numpy()}, index=adata.obs_names)
    if args.state_key:
        if args.state_key not in adata.obs:
            raise SystemExit(f"AnnData is missing obs['{args.state_key}'].")
        table["state"] = adata.obs[args.state_key].astype(str).to_numpy()
        summary = table.groupby("state")["pseudotime"].agg(["count", "mean", "median", "min", "max"]).sort_values("median")
        summary.to_csv(f"{args.out_prefix}_state_pseudotime_summary.tsv", sep="\t")
    table.to_csv(f"{args.out_prefix}_pseudotime.tsv", sep="\t")
    print(f"Wrote spatial pseudotime outputs with prefix {args.out_prefix}")


if __name__ == "__main__":
    main()
