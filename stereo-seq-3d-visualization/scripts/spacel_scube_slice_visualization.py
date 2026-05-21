#!/usr/bin/env python3
"""SPACEL/Scube-backed black-background single-slice visualization.

Template provenance:
SPACEL, DOI 10.1038/s41467-023-43220-3,
https://github.com/QuKunLab/SPACEL, original file `SPACEL/Scube/plot.py`
functions `plot_single_slice` and `plot_stacked_slices`.

This wrapper intentionally requires SPACEL/Scube. It should be run in:
`conda run -n stereo-skills-py-spacel3d python ...`
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw


DEFAULT_PALETTE = [
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#B07AA1",
    "#76B7B2",
    "#EDC948",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#1F77B4",
    "#D62728",
    "#2CA02C",
    "#9467BD",
    "#8C564B",
    "#17BECF",
]


def require_spacel_plot():
    try:
        import SPACEL
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Missing `SPACEL`. Create and use the dedicated environment:\n"
            "  bash stereo-seq-3d-visualization/envs/install-python-spacel3d.sh\n"
            "  conda run -n stereo-skills-py-spacel3d python "
            "stereo-seq-3d-visualization/scripts/spacel_scube_slice_visualization.py ...\n"
            "Blocked step: SPACEL/Scube single-slice visualization."
        ) from exc
    plot_path = Path(SPACEL.__file__).resolve().parent / "Scube" / "plot.py"
    if not plot_path.exists():
        raise SystemExit(f"Installed SPACEL is missing Scube plot module: {plot_path}")
    spec = importlib.util.spec_from_file_location("spacel_scube_plot", plot_path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load SPACEL Scube plot module from {plot_path}")
    scube_plot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scube_plot)
    return scube_plot


def load_palette(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    tab = pd.read_csv(path, sep=None, engine="python")
    if not {"label", "color"}.issubset(tab.columns):
        raise SystemExit("Palette table must contain columns: label, color")
    return dict(zip(tab["label"].astype(str), tab["color"].astype(str)))


def build_adata_from_table(df: pd.DataFrame, args: argparse.Namespace):
    try:
        import anndata as ad
    except Exception as exc:  # pragma: no cover - environment guard
        raise SystemExit("Missing `anndata` in the SPACEL 3D environment.") from exc
    missing = {args.x_col, args.y_col, args.label_col}.difference(df.columns)
    if missing:
        raise SystemExit(f"Input is missing columns: {', '.join(sorted(missing))}")
    xy = df[[args.x_col, args.y_col]].apply(pd.to_numeric, errors="coerce").to_numpy()
    labels = df[args.label_col].astype(str)
    keep = np.isfinite(xy).all(axis=1) & labels.notna().to_numpy()
    xy = xy[keep]
    labels = labels.to_numpy()[keep]
    if args.invert_y:
        xy[:, 1] = -xy[:, 1]
    adata = ad.AnnData(X=np.zeros((xy.shape[0], 1), dtype=np.float32))
    adata.obsm[args.spatial_key] = xy
    adata.obs[args.label_col] = pd.Categorical(labels)
    categories = list(adata.obs[args.label_col].cat.categories)
    palette = load_palette(args.palette)
    for idx, label in enumerate(categories):
        palette.setdefault(label, DEFAULT_PALETTE[idx % len(DEFAULT_PALETTE)])
    adata.uns[f"{args.label_col}_colors"] = [palette[label] for label in categories]
    return adata, palette


def read_inputs(args: argparse.Namespace) -> list[tuple[str, Any, dict[str, str]]]:
    scube_plot = require_spacel_plot()
    _ = scube_plot
    inputs: list[tuple[str, Any, dict[str, str]]] = []
    for input_path in args.input:
        path = Path(input_path)
        if path.suffix.lower() == ".h5ad":
            try:
                import anndata as ad
            except Exception as exc:  # pragma: no cover - environment guard
                raise SystemExit("Missing `anndata` for .h5ad input in the SPACEL 3D environment.") from exc
            adata = ad.read_h5ad(path)
            if args.spatial_key not in adata.obsm:
                raise SystemExit(f"{path} is missing obsm['{args.spatial_key}'].")
            if args.label_col not in adata.obs:
                raise SystemExit(f"{path} is missing obs['{args.label_col}'].")
            adata.obs[args.label_col] = pd.Categorical(adata.obs[args.label_col].astype(str))
            categories = list(adata.obs[args.label_col].cat.categories)
            palette = load_palette(args.palette)
            if f"{args.label_col}_colors" not in adata.uns:
                for idx, label in enumerate(categories):
                    palette.setdefault(label, DEFAULT_PALETTE[idx % len(DEFAULT_PALETTE)])
                adata.uns[f"{args.label_col}_colors"] = [palette[label] for label in categories]
            else:
                palette.update(dict(zip(categories, adata.uns[f"{args.label_col}_colors"])))
            if args.invert_y:
                adata = adata.copy()
                adata.obsm[args.spatial_key] = np.asarray(adata.obsm[args.spatial_key]).copy()
                adata.obsm[args.spatial_key][:, 1] = -adata.obsm[args.spatial_key][:, 1]
            inputs.append((path.stem, adata, palette))
            continue
        df = pd.read_csv(path, sep=None, engine="python")
        if args.slice_col and args.slice_col in df.columns:
            for slice_name, sub in df.groupby(args.slice_col, sort=True):
                adata, palette = build_adata_from_table(sub, args)
                inputs.append((f"{path.stem}_{slice_name}", adata, palette))
        else:
            adata, palette = build_adata_from_table(df, args)
            inputs.append((path.stem, adata, palette))
    return inputs


def style_dark(fig: Any) -> None:
    fig.patch.set_facecolor("black")
    for ax in fig.axes:
        ax.set_facecolor("black")
        for spine in ax.spines.values():
            spine.set_color("black")
        ax.tick_params(colors="white")


def render_single(scube_plot: Any, name: str, adata: Any, args: argparse.Namespace) -> Path:
    fig = plt.figure(figsize=(args.width, args.height))
    scube_plot.plot_single_slice(
        adata,
        spatial_key=args.spatial_key,
        cluster_key=args.label_col,
        frameon=not args.hide_axes,
        i=1,
        j=1,
        n=1,
        s=args.point_size,
    )
    style_dark(fig)
    ax = fig.axes[0]
    ax.set_aspect("equal")
    if args.hide_axes:
        ax.set_axis_off()
    if args.title:
        ax.set_title(args.title if len(args.input) == 1 else name, color="white", fontsize=11)
    out = args.outdir / f"{args.prefix}_{name}.png"
    fig.savefig(out, dpi=args.dpi, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.02)
    if args.pdf:
        fig.savefig(out.with_suffix(".pdf"), facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return out


def make_contact_sheet(images: list[Path], out: Path, cols: int = 4) -> None:
    loaded = [Image.open(path).convert("RGB") for path in images]
    if not loaded:
        return
    width = max(img.width for img in loaded)
    height = max(img.height for img in loaded)
    rows = int(math.ceil(len(loaded) / cols))
    sheet = Image.new("RGB", (cols * width, rows * height), "black")
    draw = ImageDraw.Draw(sheet)
    for idx, (path, img) in enumerate(zip(images, loaded)):
        x = (idx % cols) * width
        y = (idx // cols) * height
        sheet.paste(img, (x, y))
        draw.text((x + 12, y + 10), path.stem, fill="white")
    sheet.save(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, help="CSV/TSV or .h5ad. May be repeated.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", default="spacel_slice")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--slice-col", default="slice")
    parser.add_argument("--label-col", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--palette", default="")
    parser.add_argument("--contact-sheet", action="store_true")
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--width", type=float, default=5.2)
    parser.add_argument("--height", type=float, default=5.2)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--pdf", action="store_true")
    parser.add_argument("--hide-axes", action="store_true", default=True)
    parser.add_argument("--show-axes", dest="hide_axes", action="store_false")
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--title", default="")
    args = parser.parse_args()

    scube_plot = require_spacel_plot()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rendered: list[Path] = []
    palettes: dict[str, dict[str, str]] = {}
    for name, adata, palette in read_inputs(args):
        rendered.append(render_single(scube_plot, name, adata, args))
        palettes[name] = palette
    contact_path = ""
    if args.contact_sheet or len(rendered) > 1:
        contact = args.outdir / f"{args.prefix}_contact_sheet.png"
        make_contact_sheet(rendered, contact)
        contact_path = str(contact)
    provenance = {
        "tool": "SPACEL.Scube.plot.plot_single_slice",
        "paper": "SPACEL: deep learning-based characterization of spatial transcriptome architectures",
        "doi": "10.1038/s41467-023-43220-3",
        "repository": "https://github.com/QuKunLab/SPACEL",
        "source_files": ["SPACEL/Scube/plot.py"],
        "inputs": args.input,
        "label_col": args.label_col,
        "spatial_key": args.spatial_key,
        "black_background": True,
        "hide_axes": args.hide_axes,
        "rendered": [str(path) for path in rendered],
        "contact_sheet": contact_path,
        "palettes": palettes,
    }
    provenance_path = args.outdir / f"{args.prefix}_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Wrote {len(rendered)} SPACEL/Scube slice panel(s) to {args.outdir}")
    if contact_path:
        print(f"Wrote contact sheet: {contact_path}")
    print(f"Wrote provenance: {provenance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
