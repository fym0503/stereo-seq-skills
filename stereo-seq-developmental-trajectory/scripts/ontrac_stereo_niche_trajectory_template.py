#!/usr/bin/env python3
"""ONTraC Stereo-seq niche trajectory command template.

Template provenance:
Paper: ONTraC characterizes spatially continuous variations of tissue
microenvironment through niche trajectory analysis
Paper DOI: 10.1186/s13059-025-03588-5
Code DOI: 10.5281/zenodo.14171604
Original code: https://github.com/gyuanlab/ONTraC_paper,
`Stereo_seq_midbrain_data/run_ONTraC/stereo_midbrain_base_run_lsf.sh` and
`stereo_midbrain_base_analysis_lsf.sh`.

Reusable success pattern: keep ONTraC training and ONTraC_analysis commands
explicit, write logs, preserve NN/GNN/NT output directories, and avoid silently
falling back when the ONTraC CLI is unavailable.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def require_cli(command: str) -> str:
    path = shutil.which(command)
    if path is None:
        raise SystemExit(
            f"Missing command `{command}`. Install or activate the ONTraC environment "
            f"before running. Blocked step: {command} niche trajectory."
        )
    return path


def run_command(command: list[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        raise SystemExit(f"Command failed with exit code {proc.returncode}. See log: {log_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta-input", required=True, help="ONTraC input CSV")
    parser.add_argument("--outdir", default="ontrac_stereo_out")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--n-cpu", type=int, default=4)
    parser.add_argument("--n-neighbors", type=int, default=50)
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--patience", type=int, default=100)
    parser.add_argument("--min-delta", type=float, default=0.001)
    parser.add_argument("--min-epochs", type=int, default=50)
    parser.add_argument("--lr", type=float, default=0.03)
    parser.add_argument("--hidden-feats", type=int, default=4)
    parser.add_argument("--n-gcn-layers", type=int, default=2)
    parser.add_argument("-k", "--k", type=int, default=6)
    parser.add_argument("--modularity-loss-weight", type=float, default=0.3)
    parser.add_argument("--purity-loss-weight", type=float, default=300)
    parser.add_argument("--regularization-loss-weight", type=float, default=0.1)
    parser.add_argument("--beta", type=float, default=0.03)
    parser.add_argument("--skip-run", action="store_true", help="Only run ONTraC_analysis on existing NN/GNN/NT dirs")
    args = parser.parse_args()

    require_cli("ONTraC")
    require_cli("ONTraC_analysis")
    outdir = Path(args.outdir)
    nn_dir = outdir / "NN"
    gnn_dir = outdir / "GNN"
    nt_dir = outdir / "NT"
    analysis_dir = outdir / "analysis"
    log_dir = outdir / "logs"
    for path in [nn_dir, gnn_dir, nt_dir, analysis_dir, log_dir]:
        path.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        run_cmd = [
            "ONTraC",
            "--meta-input",
            args.meta_input,
            "--NN-dir",
            str(nn_dir),
            "--GNN-dir",
            str(gnn_dir),
            "--NT-dir",
            str(nt_dir),
            "--n-cpu",
            str(args.n_cpu),
            "--n-neighbors",
            str(args.n_neighbors),
            "--device",
            args.device,
            "--epochs",
            str(args.epochs),
            "--batch-size",
            str(args.batch_size),
            "-s",
            str(args.seed),
            "--patience",
            str(args.patience),
            "--min-delta",
            str(args.min_delta),
            "--min-epochs",
            str(args.min_epochs),
            "--lr",
            str(args.lr),
            "--hidden-feats",
            str(args.hidden_feats),
            "--n-gcn-layers",
            str(args.n_gcn_layers),
            "-k",
            str(args.k),
            "--modularity-loss-weight",
            str(args.modularity_loss_weight),
            "--purity-loss-weight",
            str(args.purity_loss_weight),
            "--regularization-loss-weight",
            str(args.regularization_loss_weight),
            "--beta",
            str(args.beta),
        ]
        run_command(run_cmd, log_dir / "ontrac_run.log")

    analysis_cmd = [
        "ONTraC_analysis",
        "--meta-input",
        args.meta_input,
        "--NN-dir",
        str(nn_dir),
        "--GNN-dir",
        str(gnn_dir),
        "--NT-dir",
        str(nt_dir),
        "-o",
        str(analysis_dir),
        "-l",
        str(log_dir / "ontrac_analysis_internal.log"),
        "-s",
        "--suppress-cell-type-composition",
    ]
    run_command(analysis_cmd, log_dir / "ontrac_analysis.log")
    print(f"Wrote ONTraC outputs to {outdir}")


if __name__ == "__main__":
    main()
