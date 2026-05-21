#!/usr/bin/env python3
"""pySCENIC full regulon inference pipeline template.

Template provenance:
Stereo-seq corpus evidence: pySCENIC/SCENIC appears in 25 local Stereo-seq
papers, including mouse placentation (DOI 10.1038/s41421-024-00740-6) and
mouse organogenesis/STcomm (DOI 10.1038/s41467-023-40155-7).
Representative tool source: https://github.com/aertslab/pySCENIC.

Reusable success pattern: keep GRN inference, motif pruning, and AUCell as
separate auditable files; require an explicit TF list, cisTarget rankings, and
motif annotations; then reuse the plotting template for spatial regulon maps.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def require_command(command: str) -> str:
    path = shutil.which(command)
    if path is None:
        raise SystemExit(
            f"Missing command `{command}`. Install or activate a pySCENIC "
            f"environment before running. Blocked step: pySCENIC pipeline."
        )
    return path


def run(command: list[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        log.write(" ".join(command) + "\n\n")
        proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        raise SystemExit(f"Command failed with exit code {proc.returncode}. See log: {log_path}")


def export_expression_from_h5ad(h5ad: str, out_csv: Path, layer: str | None) -> None:
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing Python package `scanpy`. It is required to export expression "
            "from .h5ad for pySCENIC. Blocked step: AnnData expression export."
        ) from exc
    adata = sc.read_h5ad(h5ad)
    matrix = adata.layers[layer] if layer else adata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    df = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    df.to_csv(out_csv)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, help="Expression CSV/TSV/loom/h5ad; cells x genes for CSV/TSV")
    parser.add_argument("--tf-list", required=True)
    parser.add_argument("--ranking-db", required=True, nargs="+", help="cisTarget .feather database(s)")
    parser.add_argument("--motif-annotations", required=True)
    parser.add_argument("--outdir", default="pyscenic_pipeline_out")
    parser.add_argument("--prefix", default="stereo_pyscenic")
    parser.add_argument("--layer", default="", help="AnnData layer to export when --expression is h5ad")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-grn", action="store_true", help="Use existing adjacencies.tsv in outdir")
    args = parser.parse_args()

    require_command("pyscenic")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_dir = outdir / "logs"
    log_dir.mkdir(exist_ok=True)

    expression = Path(args.expression)
    expr_for_pyscenic = expression
    if expression.suffix.lower() == ".h5ad":
        expr_for_pyscenic = outdir / f"{args.prefix}_expression.csv"
        export_expression_from_h5ad(str(expression), expr_for_pyscenic, args.layer or None)

    adj = outdir / f"{args.prefix}_adjacencies.tsv"
    regulons = outdir / f"{args.prefix}_regulons.csv"
    auc = outdir / f"{args.prefix}_auc_mtx.csv"

    if not args.skip_grn:
        run(
            [
                "pyscenic",
                "grn",
                str(expr_for_pyscenic),
                args.tf_list,
                "-o",
                str(adj),
                "--num_workers",
                str(args.workers),
                "--seed",
                str(args.seed),
            ],
            log_dir / "01_grn.log",
        )
    elif not adj.exists():
        raise SystemExit(f"--skip-grn was set but adjacency file does not exist: {adj}")

    run(
        [
            "pyscenic",
            "ctx",
            str(adj),
            *args.ranking_db,
            "--annotations_fname",
            args.motif_annotations,
            "--expression_mtx_fname",
            str(expr_for_pyscenic),
            "--output",
            str(regulons),
            "--num_workers",
            str(args.workers),
        ],
        log_dir / "02_ctx.log",
    )
    run(
        [
            "pyscenic",
            "aucell",
            str(expr_for_pyscenic),
            str(regulons),
            "--output",
            str(auc),
            "--num_workers",
            str(args.workers),
        ],
        log_dir / "03_aucell.log",
    )

    print(f"Wrote pySCENIC outputs to {outdir}")


if __name__ == "__main__":
    main()
