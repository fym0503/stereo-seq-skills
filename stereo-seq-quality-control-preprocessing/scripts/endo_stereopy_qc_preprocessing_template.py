#!/usr/bin/env python3
"""StereoPy GEF QC and preprocessing template.

Template provenance:
Paper: Single-cell profiling of the human endometrium in polycystic ovary syndrome
Paper DOI: 10.1038/s41591-025-03592-z
Original code: https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R,
`Spatial RNA-seq analysis/00.1_Stereopy_binning.py`.

Reusable success pattern: load one Stereo-seq GEF per run, set bin size
explicitly, calculate QC before and after filtering, save a filtered h5ad
checkpoint for downstream Seurat/StereoPy analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def require_stereo():
    try:
        import stereo as st
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `stereo`/Stereopy. Activate the environment "
            "containing Stereopy before running. Blocked step: GEF QC preprocessing."
        ) from exc
    return st


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gef", required=True)
    parser.add_argument("--outdir", default="stereopy_qc_preprocess")
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--min-gene", type=int, default=20)
    parser.add_argument("--min-n-genes-by-counts", type=int, default=250)
    parser.add_argument("--max-pct-mt", type=float, default=20)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--normalize-log1p", action="store_true")
    args = parser.parse_args()

    st = require_stereo()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = st.io.read_gef(file_path=args.gef, bin_size=args.bin_size)
    data.tl.cal_qc()
    data.plt.violin(out_path=str(outdir / f"{args.sample}_qc_violin_before_filter.pdf"))
    data.plt.genes_count(out_path=str(outdir / f"{args.sample}_gene_count_before_filter.pdf"))
    data.plt.spatial_scatter(out_path=str(outdir / f"{args.sample}_spatial_qc_before_filter.pdf"))

    data.tl.filter_cells(
        min_gene=args.min_gene,
        min_n_genes_by_counts=args.min_n_genes_by_counts,
        pct_counts_mt=args.max_pct_mt,
        inplace=True,
    )
    data.tl.cal_qc()
    data.plt.violin(out_path=str(outdir / f"{args.sample}_qc_violin_after_filter.pdf"))
    data.plt.genes_count(out_path=str(outdir / f"{args.sample}_gene_count_after_filter.pdf"))
    data.plt.spatial_scatter(out_path=str(outdir / f"{args.sample}_spatial_qc_after_filter.pdf"))

    st.io.write_h5ad(data, output=str(outdir / f"{args.sample}_filtered_bin{args.bin_size}.h5ad"))
    if args.normalize_log1p:
        data.tl.raw_checkpoint()
        data.tl.normalize_total(target_sum=args.target_sum)
        data.tl.log1p()
        st.io.write_h5ad(data, output=str(outdir / f"{args.sample}_filtered_bin{args.bin_size}_log1p.h5ad"))

    print(f"Wrote StereoPy QC/preprocessing outputs to {outdir}")


if __name__ == "__main__":
    main()
