#!/usr/bin/env python3
"""Spatial cell-type/domain enrichment template.

Template provenance:
P09 AD prefrontal cortex, DOI 10.1038/s41467-024-54715-y,
Zenodo 10.5281/zenodo.14048103, original files `2_layer_annotation.ipynb`,
`5_location_Tau.py`, `case_tau_location.py`, and `control_tau_location.py`.
Human cortex single-cell resolution atlas, DOI 10.1038/s41467-025-62793-9,
https://github.com/lcy1364/Cortex-Atlas-Code, original files
`src/STEREO/2_Deconvolution_and_QC/4_spatialCellMeta.R` and
`src/STEREO/4_cellSomaProximity/somafrequentgraph_all.R`.
Avian optic tectum atlas, DOI 10.1016/j.isci.2024.109009,
https://github.com/Coleliao/Spatial_OT, original cellbin/spatial annotation
interpretation scripts.

Reusable success pattern: quantify enrichment before making spatial claims,
show observed/expected ratios as a readable heatmap, and pair the table with
clear caveats about annotation and coordinate resolution.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta", required=True, help="CSV/TSV metadata table.")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--domain-col", default="domain")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-count", type=int, default=5)
    parser.add_argument("--width", type=float, default=7.0)
    parser.add_argument("--height", type=float, default=5.4)
    args = parser.parse_args()

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    meta = read_table(args.meta)
    missing = {args.celltype_col, args.domain_col}.difference(meta.columns)
    if missing:
        raise SystemExit(f"Metadata table is missing columns: {', '.join(sorted(missing))}")

    tab = pd.crosstab(meta[args.domain_col].astype(str), meta[args.celltype_col].astype(str))
    tab = tab.loc[tab.sum(axis=1) >= args.min_count, tab.sum(axis=0) >= args.min_count]
    expected = np.outer(tab.sum(axis=1), tab.sum(axis=0)) / max(tab.to_numpy().sum(), 1)
    oe = tab.to_numpy() / np.maximum(expected, 1e-12)
    oe_df = pd.DataFrame(oe, index=tab.index, columns=tab.columns)
    frac_df = tab.div(tab.sum(axis=1), axis=0)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tab.to_csv(f"{args.out_prefix}_counts.tsv", sep="\t")
    oe_df.to_csv(f"{args.out_prefix}_observed_expected.tsv", sep="\t")
    frac_df.to_csv(f"{args.out_prefix}_domain_fractions.tsv", sep="\t")

    plot_values = np.log2(oe_df.replace(0, np.nan)).replace([np.inf, -np.inf], np.nan).fillna(0)
    vmax = max(1.0, float(np.nanpercentile(np.abs(plot_values.to_numpy()), 98)))
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    im = ax.imshow(plot_values.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(plot_values.shape[1]))
    ax.set_xticklabels(plot_values.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(plot_values.shape[0]))
    ax.set_yticklabels(plot_values.index, fontsize=9)
    ax.set_xlabel(args.celltype_col, fontsize=10)
    ax.set_ylabel(args.domain_col, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2 observed/expected", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_enrichment_heatmap.pdf", dpi=300)
    plt.close(fig)
    print(f"Wrote spatial cell-type enrichment outputs with prefix {args.out_prefix}")


if __name__ == "__main__":
    main()
