#!/usr/bin/env python3
"""Stereopy spatial QC, clustering, marker, and domain plotting template.

Template provenance:
Paper: Single-cell profiling of the human endometrium in polycystic ovary syndrome
Paper DOI: 10.1038/s41591-025-03592-z
Original code: ReproductiveEndocrinologyMetabolism/Endo.R,
  Spatial RNA-seq analysis/00.1_Stereopy_binning.py
Reusable success pattern: QC maps, spatial clustering, marker-gene spatial
scatter, and cluster scatter generated from one Stereo-seq object.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def require_stereo():
    try:
        import stereo as st
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `stereopy`/`stereo`. Install or activate the "
            "environment containing Stereopy before running this template. "
            "Blocked step: Stereopy spatial QC/domain plotting."
        ) from exc
    return st


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gef", required=True, help="Stereo-seq GEF input")
    parser.add_argument("--outdir", default="stereopy_domain_out")
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--markers", default="", help="Comma-separated marker genes to plot")
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()

    st = require_stereo()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = st.io.read_gef(args.gef, bin_size=args.bin_size)
    data.tl.cal_qc()
    data.plt.spatial_scatter(out_path=str(outdir / "qc_spatial_scatter.pdf"))

    data.tl.raw_checkpoint()
    data.tl.normalize_total(target_sum=1e4)
    data.tl.log1p()
    data.tl.highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5)
    data.tl.pca(use_highly_genes=True, n_pcs=30)
    data.tl.neighbors(pca_res_key="pca", n_pcs=30)
    data.tl.umap()
    data.tl.leiden(resolution=args.resolution, res_key="leiden")
    data.plt.cluster_scatter(res_key="leiden", out_path=str(outdir / "leiden_spatial_domains.pdf"))

    markers = [gene.strip() for gene in args.markers.split(",") if gene.strip()]
    if markers:
        data.plt.spatial_scatter_by_gene(gene_name=markers, out_path=str(outdir / "marker_spatial_maps.pdf"))
        data.tl.find_marker_genes(cluster_res_key="leiden", method="t_test")
        data.plt.marker_genes_scatter(res_key="marker_genes", genes=markers, out_path=str(outdir / "marker_dot_scatter.pdf"))

    st.io.write_h5ad(data, output=str(outdir / "stereopy_domain_result.h5ad"))
    print(f"Wrote Stereopy domain outputs to {outdir}")


if __name__ == "__main__":
    main()
