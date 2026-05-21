#!/usr/bin/env python3
"""Highlight target cells on tissue coordinates using the P09 red/grey pattern.

Source: adapted from P09 Zenodo `5_location_Tau.py`, `case_tau_location.py`,
and `control_tau_location.py`.
Paper DOI: 10.1038/s41467-024-54715-y
Code DOI: 10.5281/zenodo.14048103, license recorded by Zenodo: CC-BY-4.0.
Reused successful pattern: focal target cells in red, all other cells in grey,
manual legend, optional y-axis inversion, 300 dpi output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CSV/TSV with x, y, and target indicator/label columns.")
    parser.add_argument("--x", default="x")
    parser.add_argument("--y", default="y")
    parser.add_argument("--label", required=True, help="Column used to select highlighted cells.")
    parser.add_argument("--target", required=True, help="Exact value or substring to highlight.")
    parser.add_argument("--substring", action="store_true", help="Treat --target as substring instead of exact value.")
    parser.add_argument("--out", default="p09_target_cell_location.png")
    parser.add_argument("--target-name", default="Target")
    parser.add_argument("--size", type=float, default=5)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--invert-x", action="store_true")
    args = parser.parse_args()

    df = read_table(Path(args.input))
    missing = {args.x, args.y, args.label} - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {', '.join(sorted(missing))}")

    values = df[args.label].astype(str)
    is_target = values.str.contains(args.target, regex=False, na=False) if args.substring else values.eq(args.target)
    colors = is_target.map({True: "red", False: "grey"})

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(df[args.x], df[args.y], c=colors, s=args.size, linewidths=0)
    handles = [
        Line2D([0], [0], marker="o", color="w", label=args.target_name, markersize=8, markerfacecolor="red"),
        Line2D([0], [0], marker="o", color="w", label="Others", markersize=8, markerfacecolor="grey"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=12)
    if args.invert_y:
        ax.invert_yaxis()
    if args.invert_x:
        ax.invert_xaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)


if __name__ == "__main__":
    main()
