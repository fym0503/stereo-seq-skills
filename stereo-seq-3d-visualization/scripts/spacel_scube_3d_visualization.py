#!/usr/bin/env python3
"""SPACEL/Scube-backed black-background 3D visualization for Stereo-seq.

Template provenance:
SPACEL, DOI 10.1038/s41467-023-43220-3,
https://github.com/QuKunLab/SPACEL, original files `SPACEL/Scube/plot.py`,
`SPACEL/Scube/utils_3d.py`, and `docs/tutorials/Stereo-seq_Scube.ipynb`.

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


VIEW_PRESETS = {
    "front": (12.0, -90.0),
    "side": (10.0, 0.0),
    "top": (90.0, -90.0),
    "oblique": (25.0, -60.0),
    "poster": (18.0, -45.0),
}

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
            "stereo-seq-3d-visualization/scripts/spacel_scube_3d_visualization.py ...\n"
            "Blocked step: SPACEL/Scube 3D visualization."
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


def read_input(path: Path, args: argparse.Namespace) -> pd.DataFrame:
    if path.suffix.lower() == ".h5ad":
        try:
            import anndata as ad
        except Exception as exc:  # pragma: no cover - environment guard
            raise SystemExit("Missing `anndata` for .h5ad input in the SPACEL 3D environment.") from exc
        adata = ad.read_h5ad(path)
        if args.obsm_key not in adata.obsm:
            raise SystemExit(f"AnnData is missing obsm['{args.obsm_key}'].")
        coords = np.asarray(adata.obsm[args.obsm_key])
        if coords.shape[1] < 2:
            raise SystemExit(f"obsm['{args.obsm_key}'] must contain at least x/y columns.")
        df = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=adata.obs_names)
        if coords.shape[1] >= 3:
            df["z"] = coords[:, 2]
        for col in [args.slice_col, args.value_col, args.color_col]:
            if col and col in adata.obs:
                df[col] = adata.obs[col].to_numpy()
        if args.value_col and args.value_col not in df.columns and args.value_col in adata.var_names:
            values = adata[:, args.value_col].X
            if hasattr(values, "toarray"):
                values = values.toarray()
            df[args.value_col] = np.asarray(values).ravel()
        return df.reset_index(names="unit_id")
    return pd.read_csv(path, sep=None, engine="python")


def stable_slice_order(values: pd.Series, explicit: str) -> list[str]:
    if explicit:
        order = [item.strip() for item in explicit.split(",") if item.strip()]
        missing = set(values.astype(str).unique()).difference(order)
        return order + sorted(missing)
    unique = values.astype(str).unique().tolist()
    try:
        return sorted(unique, key=lambda item: float(item))
    except ValueError:
        return sorted(unique)


def prepare_coords(df: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, pd.DataFrame]:
    missing = {args.x_col, args.y_col}.difference(df.columns)
    if missing:
        raise SystemExit(f"Input is missing coordinate columns: {', '.join(sorted(missing))}")
    out = df.copy()
    out["__x"] = pd.to_numeric(out[args.x_col], errors="coerce")
    out["__y"] = pd.to_numeric(out[args.y_col], errors="coerce")
    if args.z_col and args.z_col in out.columns:
        out["__z"] = pd.to_numeric(out[args.z_col], errors="coerce")
    else:
        if args.slice_col not in out.columns:
            raise SystemExit(f"Input needs `{args.z_col}` or slice column `{args.slice_col}`.")
        order = stable_slice_order(out[args.slice_col], args.slice_order)
        z_by_slice = {slice_name: idx * args.slice_spacing for idx, slice_name in enumerate(order)}
        out["__slice"] = out[args.slice_col].astype(str)
        out["__z"] = out["__slice"].map(z_by_slice)
    if args.invert_y:
        out["__y"] = -out["__y"]
    out = out.dropna(subset=["__x", "__y", "__z"]).copy()
    if args.max_points and len(out) > args.max_points:
        out = out.sample(n=args.max_points, random_state=args.seed).sort_index()
    coords = out[["__x", "__y", "__z"]].to_numpy(dtype=float)
    return coords, out


def load_palette(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    tab = pd.read_csv(path, sep=None, engine="python")
    if not {"label", "color"}.issubset(tab.columns):
        raise SystemExit("Palette table must contain columns: label, color")
    return dict(zip(tab["label"].astype(str), tab["color"].astype(str)))


def make_colors(df: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray | None, np.ndarray | None, dict[str, str]]:
    if args.color_col and args.color_col in df.columns:
        colors = df[args.color_col].astype(str).to_numpy()
        return None, colors, {}
    if not args.value_col or args.value_col not in df.columns:
        if "__slice" in df.columns:
            args.value_col = "__slice"
        else:
            return np.full(len(df), 1.0), None, {}
    values = df[args.value_col]
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().all():
        return numeric.to_numpy(dtype=float), None, {}
    labels = values.astype(str).to_numpy()
    palette = load_palette(args.palette)
    unique = sorted(pd.unique(labels))
    for idx, label in enumerate(unique):
        palette.setdefault(label, DEFAULT_PALETTE[idx % len(DEFAULT_PALETTE)])
    return None, np.asarray([palette[label] for label in labels]), palette


def parse_views(view_args: list[str]) -> list[tuple[str, float, float]]:
    if not view_args:
        view_args = ["poster"]
    views: list[tuple[str, float, float]] = []
    for item in view_args:
        if item in VIEW_PRESETS:
            elev, azim = VIEW_PRESETS[item]
            views.append((item, elev, azim))
            continue
        parts = [part.strip() for part in item.split(",")]
        if len(parts) != 2:
            raise SystemExit(f"View must be a preset or elev,azim pair: {item}")
        elev, azim = float(parts[0]), float(parts[1])
        views.append((f"elev{elev:g}_azim{azim:g}", elev, azim))
    return views


def style_spacel_dark(fig: Any) -> None:
    fig.patch.set_facecolor("black")
    for ax in fig.axes:
        ax.set_facecolor("black")
        ax.grid(False)
        ax.tick_params(colors="white")
        for axis in [getattr(ax, "xaxis", None), getattr(ax, "yaxis", None), getattr(ax, "zaxis", None)]:
            if axis is not None:
                axis.label.set_color("white")
                try:
                    axis.pane.set_facecolor((0, 0, 0, 1))
                    axis.pane.set_edgecolor((0, 0, 0, 1))
                except Exception:
                    pass


def render_view(
    scube_plot: Any,
    coords: np.ndarray,
    val: np.ndarray | None,
    color: np.ndarray | None,
    out: Path,
    elev: float,
    azim: float,
    args: argparse.Namespace,
) -> None:
    fig = scube_plot.plot_3d(
        coords,
        val=val,
        color=color,
        figsize=(args.width, args.height),
        return_fig=True,
        elev=elev,
        azim=azim,
        frameon=not args.hide_axes,
        save_path=None,
        show=False,
        s=args.point_size,
        alpha=args.alpha,
        cmap=args.cmap,
        linewidths=0,
        depthshade=args.depthshade,
    )
    style_spacel_dark(fig)
    if args.hide_axes:
        for ax in fig.axes:
            ax.set_axis_off()
    fig.savefig(out, dpi=args.dpi, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.02)
    if args.pdf:
        fig.savefig(out.with_suffix(".pdf"), facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def make_contact_sheet(images: list[Path], out: Path, cols: int = 3) -> None:
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


def make_rotation_gif(
    scube_plot: Any,
    coords: np.ndarray,
    val: np.ndarray | None,
    color: np.ndarray | None,
    outdir: Path,
    args: argparse.Namespace,
) -> Path:
    frame_dir = outdir / f"{args.prefix}_rotation_frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    frames: list[Path] = []
    for idx, azim in enumerate(np.linspace(args.gif_start, args.gif_stop, args.gif_frames, endpoint=False)):
        frame = frame_dir / f"frame_{idx:03d}.png"
        render_view(scube_plot, coords, val, color, frame, args.gif_elev, float(azim), args)
        frames.append(frame)
    gif_path = outdir / f"{args.prefix}_rotation.gif"
    import imageio.v2 as imageio

    images = [imageio.imread(frame) for frame in frames]
    imageio.mimsave(gif_path, images, duration=args.gif_duration)
    return gif_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="CSV/TSV or .h5ad with spatial coordinates.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", default="spacel_3d")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--z-col", default="z")
    parser.add_argument("--slice-col", default="slice")
    parser.add_argument("--slice-order", default="", help="Comma-separated slice order.")
    parser.add_argument("--slice-spacing", type=float, default=1.0)
    parser.add_argument("--value-col", default="", help="Categorical or continuous value column.")
    parser.add_argument("--color-col", default="", help="Optional precomputed color column.")
    parser.add_argument("--palette", default="", help="Optional label,color table.")
    parser.add_argument("--obsm-key", default="spatial")
    parser.add_argument("--view", action="append", default=[], help="Preset or elev,azim. May be repeated.")
    parser.add_argument("--contact-sheet", action="store_true")
    parser.add_argument("--gif", action="store_true")
    parser.add_argument("--gif-frames", type=int, default=36)
    parser.add_argument("--gif-start", type=float, default=0.0)
    parser.add_argument("--gif-stop", type=float, default=360.0)
    parser.add_argument("--gif-elev", type=float, default=18.0)
    parser.add_argument("--gif-duration", type=float, default=0.08)
    parser.add_argument("--point-size", type=float, default=3.0)
    parser.add_argument("--alpha", type=float, default=0.95)
    parser.add_argument("--cmap", default="Spectral_r")
    parser.add_argument("--width", type=float, default=8.0)
    parser.add_argument("--height", type=float, default=8.0)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--pdf", action="store_true")
    parser.add_argument("--hide-axes", action="store_true", default=True)
    parser.add_argument("--show-axes", dest="hide_axes", action="store_false")
    parser.add_argument("--invert-y", action="store_true")
    parser.add_argument("--max-points", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--depthshade", action="store_true")
    args = parser.parse_args()

    scube_plot = require_spacel_plot()
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = read_input(args.input, args)
    coords, prepared = prepare_coords(df, args)
    val, color, palette = make_colors(prepared, args)

    snapshot_paths: list[Path] = []
    for name, elev, azim in parse_views(args.view):
        out = args.outdir / f"{args.prefix}_{name}.png"
        render_view(scube_plot, coords, val, color, out, elev, azim, args)
        snapshot_paths.append(out)

    contact_path = None
    if args.contact_sheet or len(snapshot_paths) > 1:
        contact_path = args.outdir / f"{args.prefix}_contact_sheet.png"
        make_contact_sheet(snapshot_paths, contact_path)

    gif_path = None
    if args.gif:
        gif_path = make_rotation_gif(scube_plot, coords, val, color, args.outdir, args)

    provenance = {
        "tool": "SPACEL.Scube.plot.plot_3d",
        "paper": "SPACEL: deep learning-based characterization of spatial transcriptome architectures",
        "doi": "10.1038/s41467-023-43220-3",
        "repository": "https://github.com/QuKunLab/SPACEL",
        "source_files": ["SPACEL/Scube/plot.py", "SPACEL/Scube/utils_3d.py", "docs/tutorials/Stereo-seq_Scube.ipynb"],
        "input": str(args.input),
        "n_points_rendered": int(len(prepared)),
        "x_col": args.x_col,
        "y_col": args.y_col,
        "z_col": args.z_col if args.z_col in prepared.columns else "",
        "slice_col": args.slice_col if args.slice_col in prepared.columns else "",
        "slice_spacing": args.slice_spacing,
        "value_col": args.value_col,
        "color_col": args.color_col,
        "views": [{"name": n, "elev": e, "azim": a} for n, e, a in parse_views(args.view)],
        "black_background": True,
        "hide_axes": args.hide_axes,
        "palette": palette,
        "snapshots": [str(path) for path in snapshot_paths],
        "contact_sheet": str(contact_path) if contact_path else "",
        "rotation_gif": str(gif_path) if gif_path else "",
    }
    provenance_path = args.outdir / f"{args.prefix}_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Wrote {len(snapshot_paths)} SPACEL/Scube 3D snapshot(s) to {args.outdir}")
    if contact_path:
        print(f"Wrote contact sheet: {contact_path}")
    if gif_path:
        print(f"Wrote rotation GIF: {gif_path}")
    print(f"Wrote provenance: {provenance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
