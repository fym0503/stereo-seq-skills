#!/usr/bin/env python3
"""Draw ordered cortical-layer/domain maps using the P09 layer color scheme.

Source: adapted from P09 Zenodo `2_layer_annotation.ipynb`.
Paper DOI: 10.1038/s41467-024-54715-y
Code DOI: 10.5281/zenodo.14048103, license recorded by Zenodo: CC-BY-4.0.
Reused successful pattern: fixed cortical-layer palette, ordered plotting,
axis-free equal-aspect spatial map, 300 dpi output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


P09_LAYER_COLORS = {
    "L1": "#FF7F00",
    "L2/3": "#984EA3",
    "L4": "#4DAF4A",
    "L5": "#377EB8",
    "L6": "#E41A1C",
    "WM": "#A65628",
}
P09_LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CSV/TSV with x, y, and layer/domain columns.")
    parser.add_argument("--x", default="x")
    parser.add_argument("--y", default="y")
    parser.add_argument("--label", default="annotation")
    parser.add_argument("--out", default="p09_layer_spatial_map.png")
    parser.add_argument("--size", type=float, default=4)
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--invert-x", action="store_true")
    parser.add_argument("--title", default="")
    args = parser.parse_args()

    df = read_table(Path(args.input))
    missing = {args.x, args.y, args.label} - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {', '.join(sorted(missing))}")

    fig, ax = plt.subplots(figsize=(8, 8))
    labels = [label for label in P09_LAYER_ORDER if label in set(df[args.label])]
    labels += [label for label in sorted(set(df[args.label])) if label not in P09_LAYER_COLORS]

    fallback = plt.get_cmap("tab20")
    for index, label in enumerate(labels):
        sub = df[df[args.label] == label]
        color = P09_LAYER_COLORS.get(label, fallback(index % 20))
        ax.scatter(sub[args.x], sub[args.y], c=[color], s=args.size, linewidths=0, label=label)

    if args.invert_y:
        ax.invert_yaxis()
    if args.invert_x:
        ax.invert_xaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    if args.title:
        ax.set_title(args.title)
    ax.legend(frameon=False, markerscale=3, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)


if __name__ == "__main__":
    main()
