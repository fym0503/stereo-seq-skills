#!/usr/bin/env python3
"""QC template for Stereo-seq cellbin segmentation masks and coordinates.

Template provenance:
STCellbin, DOI 10.46471/gigabyte.110, https://github.com/STOmics/STCellbin,
original files `STCellbin.py`, `src/cellbin/modules/cell_segmentation.py`,
and `src/cellbin/contrib/cell_mask.py`.
BIDCell/CellSPA, DOI 10.1038/s41467-023-44560-w,
https://github.com/SydneyBioX/BIDCell and https://github.com/SydneyBioX/CellSPA,
original files `bidcell/example_params/stereoseq.yaml`,
`bidcell/processing/nuclei_segmentation.py`, and `R/visualisation.R`.
Ascidian endostyle Stereo-seq, DOI 10.1126/sciadv.adi9035,
https://github.com/lskfs/ascidian-endostyle, original files
`01.stereo-seq_cellseg/cellsegment.py`, `centroid.py`, and `gem_mask.py`.

Reusable success pattern: separate mask geometry QC from expression aggregation,
plot equal-aspect spatial/image coordinates, use readable Arial-compatible labels,
and report missing segmentation dependencies instead of silently switching methods.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def read_mask(path: Path) -> np.ndarray:
    try:
        import imageio.v3 as iio
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `imageio`. Install it before reading mask images. "
            "Blocked step: cellbin mask QC."
        ) from exc
    mask = iio.imread(path)
    if mask.ndim == 3:
        mask = mask[..., 0]
    return np.asarray(mask)


def mask_summary(mask: np.ndarray) -> pd.DataFrame:
    labels, counts = np.unique(mask, return_counts=True)
    keep = labels != 0
    labels = labels[keep]
    counts = counts[keep]
    return pd.DataFrame({"cell_id": labels.astype(str), "area_px": counts.astype(int)})


def centroid_table(mask: np.ndarray) -> pd.DataFrame:
    labels = np.unique(mask)
    labels = labels[labels != 0]
    records = []
    for label in labels:
        ys, xs = np.where(mask == label)
        if len(xs) == 0:
            continue
        records.append(
            {
                "cell_id": str(label),
                "x": float(xs.mean()),
                "y": float(ys.mean()),
                "area_px": int(len(xs)),
            }
        )
    return pd.DataFrame(records)


def plot_hist(values: pd.Series, xlabel: str, out: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.hist(values.dropna(), bins=50, color="#4E79A7", edgecolor="white", linewidth=0.3)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Cell count", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def plot_centroids(df: pd.DataFrame, out: Path, title: str, point_size: float) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.scatter(df["x"], df["y"], s=point_size, c="#4E79A7", alpha=0.9, linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def plot_expression_overlay(
    cell_df: pd.DataFrame,
    expr_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    out: Path,
    point_size: float,
) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.scatter(
        expr_df[x_col],
        expr_df[y_col],
        s=point_size,
        c="#BAB0AC",
        alpha=0.35,
        linewidths=0,
        label="expression",
    )
    ax.scatter(
        cell_df["x"],
        cell_df["y"],
        s=max(point_size * 2, 1.0),
        c="#E15759",
        alpha=0.9,
        linewidths=0,
        label="cell centroids",
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.set_title("Expression coordinates vs cell centroids", fontsize=11)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mask", help="Label mask image where 0 is background and each nonzero value is a cell id.")
    parser.add_argument("--cell-table", help="Optional table with cell_id and x/y or area columns.")
    parser.add_argument("--expression-table", help="Optional expression/bin coordinate table for overlay QC.")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--count-col", default=None, help="Optional count column in cell table.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="cellbin_qc")
    parser.add_argument("--point-size", type=float, default=1.0)
    args = parser.parse_args()

    if not args.mask and not args.cell_table:
        raise SystemExit("Provide --mask and/or --cell-table for cellbin segmentation QC.")

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tables = []
    centroids = None
    if args.mask:
        mask = read_mask(Path(args.mask))
        summary = mask_summary(mask)
        centroids = centroid_table(mask)
        summary_path = outdir / f"{args.prefix}_mask_summary.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)
        centroids.to_csv(outdir / f"{args.prefix}_mask_centroids.tsv", sep="\t", index=False)
        tables.append(summary)
        if not summary.empty:
            plot_hist(summary["area_px"], "Mask area (pixels)", outdir / f"{args.prefix}_mask_area.pdf", "Cell mask area")
        if centroids is not None and not centroids.empty:
            plot_centroids(centroids, outdir / f"{args.prefix}_centroids.pdf", "Cell mask centroids", args.point_size)

    cell_df = None
    if args.cell_table:
        cell_df = read_table(Path(args.cell_table))
        if args.x_col in cell_df.columns and args.y_col in cell_df.columns:
            centroids = cell_df.rename(columns={args.x_col: "x", args.y_col: "y"})
            plot_centroids(centroids, outdir / f"{args.prefix}_cell_table_coordinates.pdf", "Cell table coordinates", args.point_size)
        area_candidates = [col for col in ["area_px", "area", "cell_area"] if col in cell_df.columns]
        if area_candidates:
            plot_hist(pd.to_numeric(cell_df[area_candidates[0]], errors="coerce"), area_candidates[0], outdir / f"{args.prefix}_cell_area.pdf", "Cell area")
        if args.count_col and args.count_col in cell_df.columns:
            plot_hist(pd.to_numeric(cell_df[args.count_col], errors="coerce"), args.count_col, outdir / f"{args.prefix}_cell_counts.pdf", "Counts per cell")
        cell_df.to_csv(outdir / f"{args.prefix}_cell_table_copy.tsv", sep="\t", index=False)

    if args.expression_table and centroids is not None and not centroids.empty:
        expr_df = read_table(Path(args.expression_table))
        missing = {args.x_col, args.y_col}.difference(expr_df.columns)
        if missing:
            raise SystemExit(f"Expression table is missing coordinate columns: {', '.join(sorted(missing))}")
        plot_expression_overlay(
            centroids,
            expr_df,
            args.x_col,
            args.y_col,
            outdir / f"{args.prefix}_expression_overlay.pdf",
            args.point_size,
        )

    report = outdir / f"{args.prefix}_report.txt"
    with report.open("w", encoding="utf-8") as handle:
        handle.write("Cellbin segmentation QC report\n")
        handle.write(
            "Template sources: STCellbin DOI 10.46471/gigabyte.110; "
            "BIDCell DOI 10.1038/s41467-023-44560-w; "
            "ascidian endostyle DOI 10.1126/sciadv.adi9035\n"
        )
        if args.mask:
            handle.write(f"Mask: {args.mask}\n")
            detected = len(tables[0]) if tables else 0
            handle.write(f"Detected non-background cells: {detected}\n")
        if args.cell_table:
            handle.write(f"Cell table: {args.cell_table}\n")
        if args.expression_table:
            handle.write(f"Expression coordinate table: {args.expression_table}\n")
    print(f"Wrote cellbin QC outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
