---
name: stereo-seq-3d-visualization
description: Use when Stereo-seq or STOmics data needs black-background 3D spatial transcriptomics visualization, stacked slice snapshots, arbitrary 3D viewpoint snapshots, multi-angle 3D contact sheets, rotation GIFs, or single-slice dark spatial maps from registered multi-slice coordinates, especially using SPACEL/Scube, SLAT, STAligner, GEM3D, or Spateo-derived 3D workflows.
---

# Stereo-seq 3D Visualization

Use this skill after slices have 2D coordinates plus slice order, or after registration/alignment has produced 3D coordinates. This skill is for visual output, not for doing the registration itself.

## Use This For

- Black-background 3D ST figures like stacked slice point clouds, 3D atlas/domain maps, gene/module snapshots, and contact sheets.
- Arbitrary 3D camera/viewpoint snapshots and optional rotation GIFs.
- Single-slice black-background spatial maps that match 3D figure styling.
- Turning SPACEL/Scube, SLAT, STAligner, GEM3D, or Spateo-aligned coordinates into publication-ready 3D visualization panels.

## Default Requirements

- Use the SPACEL/Scube-backed scripts in `scripts/` first. Do not replace them with a low-dependency matplotlib-only implementation unless the user explicitly asks.
- Run in the dedicated environment: `conda run -n stereo-skills-py-spacel3d python ...`. Create it with this skill's bundled `envs/install-python-spacel3d.sh`. This environment installs SPACEL source files and pinned plotting/h5ad dependencies for visualization; it is not intended as the full SPACEL training/GPR environment.
- If `SPACEL` or snapshot/GIF dependencies are missing, stop the blocked step and report the missing package and environment file. Do not silently install packages into `stereo-skills-py`.
- Use this skill's bundled `envs/install-python-spacel3d-mesh.sh` only when a task specifically needs Open3D mesh or surface utilities from `SPACEL/Scube/utils_3d.py`.
- Read [source_code.md](references/source_code.md) before claiming provenance or adapting external 3D code.
- Keep slice spacing, coordinate orientation, camera angles, point size, downsampling, and color palette explicit in the command or report.
- Use black backgrounds for both figure and axes. Hide axes for presentation snapshots unless coordinate axes are part of the scientific claim.
- For single-slice maps, preserve equal aspect ratio and do not stretch tissue geometry.

## Workflow

1. Identify input type:
   - CSV/TSV with `x`, `y`, `slice`, and optional `z`, label/value, or color columns.
   - `.h5ad` with `obsm['spatial']` or 3D coordinate key plus `obs` label/value columns.
   - Per-slice `.h5ad` files when using SPACEL-style stacked slice panels.
2. Decide visual output:
   - static 3D snapshots;
   - multi-angle contact sheet;
   - rotation GIF;
   - black single-slice spatial map;
   - optional mesh/surface work from Scube/Open3D utilities when a dense 3D point cloud exists.
3. Use `scripts/spacel_scube_3d_visualization.py` for 3D snapshots/contact sheets/GIFs.
4. Use `scripts/spacel_scube_slice_visualization.py` for black-background single-slice maps and per-slice contact sheets.
5. Inspect output images for non-empty data, correct orientation, no legend/data overlap, and consistent colors across snapshots and slices.
6. Report provenance: SPACEL/Scube paper DOI and repository, plus any upstream alignment source such as SLAT, STAligner, GEM3D, or Spateo if those coordinates were used.

## Reusable Article Code

- `scripts/spacel_scube_3d_visualization.py`: wraps `SPACEL.Scube.plot.plot_3d` and applies black-background 3D snapshot/contact-sheet/GIF styling.
- `scripts/spacel_scube_slice_visualization.py`: prepares AnnData/table inputs for `SPACEL.Scube.plot.plot_single_slice`-style dark single-slice panels.

## Output Expectations

- PNG/PDF static snapshots for requested camera angles.
- Optional multi-angle contact sheet PNG/PDF.
- Optional rotation GIF.
- Optional single-slice black-background PNG/PDF panels.
- A small parameters/provenance TSV or JSON recording input columns, slice spacing, camera angles, downsampling, point size, palette, tool versions, and paper/code provenance.
