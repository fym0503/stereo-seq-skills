#!/usr/bin/env python3
"""Spatial ligand-receptor neighbor filtering template.

Template provenance:
STcomm mouse organogenesis, DOI 10.1038/s41467-023-40155-7,
https://github.com/gpenglab/STcomm, original files `R/cellColocation.R` and
`R/st_comm.R`.
SPIDER, DOI 10.1038/s41467-025-62988-0,
https://github.com/deepomicslab/SPIDER and SPIDER-paper notebooks.
P09 AD prefrontal cortex, DOI 10.1038/s41467-024-54715-y,
Zenodo 10.5281/zenodo.14048103, original layer/concentric CCI notebooks.

Reusable success pattern: require spatial adjacency before interpreting LR
pairs, export sender-receiver support tables, and keep LR results as hypotheses
unless expression, proximity, and biological context all agree.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def require_neighbors():
    try:
        from sklearn.neighbors import NearestNeighbors
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scikit-learn`. Install it before spatial "
            "neighbor LR filtering. Blocked step: neighbor graph construction."
        ) from exc
    return NearestNeighbors


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def expression_lookup(expr: pd.DataFrame, genes: set[str], id_col: str) -> pd.DataFrame:
    needed = [id_col] + [gene for gene in genes if gene in expr.columns]
    if len(needed) == 1:
        raise SystemExit("None of the LR genes were found as expression columns.")
    values = expr[needed].copy()
    values = values.set_index(id_col)
    return values.apply(pd.to_numeric, errors="coerce").fillna(0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta", required=True, help="CSV/TSV with id, x, y, and label columns.")
    parser.add_argument("--expr", required=True, help="CSV/TSV wide expression table keyed by id.")
    parser.add_argument("--lr", required=True, help="CSV/TSV ligand-receptor table.")
    parser.add_argument("--id-col", default="cell_id")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--label-col", default="celltype")
    parser.add_argument("--ligand-col", default="ligand")
    parser.add_argument("--receptor-col", default="receptor")
    parser.add_argument("--radius", type=float, default=100.0)
    parser.add_argument("--min-expression", type=float, default=0.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    NearestNeighbors = require_neighbors()
    meta = read_table(args.meta)
    expr = read_table(args.expr)
    lr = read_table(args.lr)
    for col in [args.id_col, args.x_col, args.y_col, args.label_col]:
        if col not in meta.columns:
            raise SystemExit(f"Metadata table missing column: {col}")
    for col in [args.id_col]:
        if col not in expr.columns:
            raise SystemExit(f"Expression table missing column: {col}")
    for col in [args.ligand_col, args.receptor_col]:
        if col not in lr.columns:
            raise SystemExit(f"LR table missing column: {col}")

    meta = meta.dropna(subset=[args.id_col, args.x_col, args.y_col, args.label_col]).copy()
    meta[args.id_col] = meta[args.id_col].astype(str)
    expr[args.id_col] = expr[args.id_col].astype(str)
    lr[args.ligand_col] = lr[args.ligand_col].astype(str)
    lr[args.receptor_col] = lr[args.receptor_col].astype(str)
    genes = set(lr[args.ligand_col]).union(set(lr[args.receptor_col]))
    expr_mat = expression_lookup(expr, genes, args.id_col)

    common = meta[args.id_col].isin(expr_mat.index)
    meta = meta[common].reset_index(drop=True)
    if meta.empty:
        raise SystemExit("No shared ids between metadata and expression table.")

    xy = meta[[args.x_col, args.y_col]].to_numpy(float)
    nn = NearestNeighbors(radius=args.radius)
    nn.fit(xy)
    neighbors = nn.radius_neighbors(xy, return_distance=False)
    id_to_label = dict(zip(meta[args.id_col], meta[args.label_col].astype(str)))
    ids = meta[args.id_col].tolist()

    records: list[dict[str, object]] = []
    for _, pair in lr.iterrows():
        ligand = pair[args.ligand_col]
        receptor = pair[args.receptor_col]
        if ligand not in expr_mat.columns or receptor not in expr_mat.columns:
            continue
        ligand_expr = expr_mat[ligand]
        receptor_expr = expr_mat[receptor]
        support: dict[tuple[str, str], list[float]] = {}
        for i, nb in enumerate(neighbors):
            sender_id = ids[i]
            if ligand_expr.get(sender_id, 0) <= args.min_expression:
                continue
            for j in nb:
                if i == j:
                    continue
                receiver_id = ids[j]
                if receptor_expr.get(receiver_id, 0) <= args.min_expression:
                    continue
                key = (id_to_label[sender_id], id_to_label[receiver_id])
                support.setdefault(key, []).append(float(ligand_expr[sender_id] * receptor_expr[receiver_id]))
        for (sender, receiver), scores in support.items():
            records.append(
                {
                    "ligand": ligand,
                    "receptor": receptor,
                    "sender": sender,
                    "receiver": receiver,
                    "neighbor_pairs": len(scores),
                    "mean_lr_product": float(np.mean(scores)),
                    "radius": args.radius,
                }
            )

    columns = [
        "ligand",
        "receptor",
        "sender",
        "receiver",
        "neighbor_pairs",
        "mean_lr_product",
        "radius",
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(records, columns=columns)
    if not result.empty:
        result = result.sort_values(["neighbor_pairs", "mean_lr_product"], ascending=False)
    result.to_csv(out, sep="\t", index=False)
    print(f"Wrote spatial LR neighbor support table to {out}")


if __name__ == "__main__":
    main()
