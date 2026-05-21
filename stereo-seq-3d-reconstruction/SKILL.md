---
name: stereo-seq-3d-reconstruction
description: Use when Stereo-seq or STOmics serial sections need multi-slice registration, 3D coordinate reconstruction, slice-to-atlas alignment, shape/contour matching, or 3D-ready visualization and mapping.
---

# Stereo-seq 3D Reconstruction

## Use This For

- Aligning serial Stereo-seq slices or histology masks into a shared coordinate frame.
- Mapping spatial coordinates through thin-plate-spline or other slice-registration transforms.
- Preparing registered coordinates for 3D reconstruction, atlas projection, or cross-section comparison.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing new registration or plotting code.
- Read [source_code.md](references/source_code.md) first for curated templates; if no curated entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files.
- Inspect local Python/R environments first. Prefer `conda run -n stereo-skills-py python ...` for Python scripts and `conda run -n stereo-skills-r Rscript ...` for R scripts. spatiAlign is configured in `stereo-skills-py`; it can also be isolated with `envs/environment-python-spatialign.yml`. If a dependency is missing, stop that step and report the missing package and blocked script.
- Use equal-aspect spatial plots, readable Arial text, non-overlapping legends, and export PDF/300 dpi diagnostic figures.
- In the final response, state the reused paper, DOI, code repository/source file, and what was changed for the current dataset.

## Workflow

1. Identify section order, coordinate units, tissue masks, anchors, and whether registration is pairwise or to a shared atlas.
2. Read [source_code.md](references/source_code.md), then adapt the closest local script:
   - `scripts/zebrafish_tps_slice_registration_template.py` for mask/anchor-based TPS coordinate mapping and QC plots.
   - `scripts/spatialign_multislice_alignment_template.py` for spatiAlign-style multi-slice spatial alignment when that package is installed.
3. Export transformed coordinates and a registration QC plot for every slice pair or atlas mapping.
4. Validate registration with landmarks, tissue contours, known anatomical layers, and non-overlap of impossible tissue regions.
5. Keep biological interpretation separate from registration quality.

## Output Expectations

- Registration method and required inputs.
- Transformed coordinate tables.
- QC plots with anchors/contours and equal aspect ratio.
- Caveats for scale, orientation, missing tissue, and manual anchors.
- Reused article code source, paper DOI, repository URL, original file name, and dataset-specific edits.
