---
name: stereo-seq-publication-plotting
description: Use when Stereo-seq or STOmics analysis needs paper-quality figures, spatial maps, dot plots, heatmaps, marker panels, module-score maps, legends, palettes, Arial typography, non-overlapping labels, or reuse of plotting code from public Stereo-seq article repositories.
---

# Stereo-seq Publication Plotting

## Use This For

- Turning Stereo-seq coordinates, annotations, gene scores, cell-type abundances, domains, or interaction scores into manuscript-ready figures.
- Reusing high-quality plotting templates from public Stereo-seq article code instead of writing figure code from scratch.
- Fixing figure issues such as tiny fonts, legends covering data, unstable palettes, non-equal spatial aspect, or raster-only output.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing new plotting code or searching external repositories.
- Read [source_code.md](references/source_code.md) when deciding which template to adapt. Let the current dataset, plot type, tissue geometry, and task similarity guide the choice; do not hard-code tissue-to-tool rules.
- If no curated plotting entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked figure scripts and reusable files before external search.
- Before running Python or R, inspect the local environment. Prefer `conda run -n stereo-skills-py python ...` for Python plotting scripts and `conda run -n stereo-skills-r Rscript ...` for R plotting scripts. If a required package is missing, stop that plotting step and tell the user which package is missing and which figure is blocked.
- Use Arial where available. For manuscript panels, keep axis text, legend text, and labels at least 9 pt by default; use 10-11 pt for legend titles, panel titles, and important labels.
- Export vector PDF plus 300 dpi raster when useful.
- Keep legends outside the data area where possible, and never let legends or colorbars cover the tissue/map content.
- In the final response, state the reused paper, DOI, code repository or code DOI, original file, and dataset-specific edits.

## Workflow

1. Identify the figure contract: input file/object, coordinate columns, grouping or value column, plot type, expected panel size, and output filenames.
2. Read [source_code.md](references/source_code.md) and choose the closest template by paper-code similarity.
3. Adapt only dataset-specific parts: file paths, column names, marker lists, labels, palette entries, bin/cell size, and output names.
4. Check that the result has equal spatial aspect, readable text, and no legend or colorbar overlap.
5. Report provenance and modifications.

## Reusable Article Code

- `scripts/stereo_spatial_panel_template.py`: categorical/continuous spatial-map template adapted from P09 layer maps, Endo.R spatial scatter patterns, and GF/SPF cecum marker spatial maps.
- `scripts/stereo_dotplot_template.R`: marker/program dotplot template adapted from Endo.R marker dotplots and GF/SPF cecum function dotplots.
- `scripts/stereo_marker_spatial_grid_template.py`: marker spatial grid template adapted from P09, Endo.R, and avian optic tectum marker-map figures.
- `scripts/stereo_marker_heatmap_template.R`: marker/program heatmap template adapted from SpaSEG, human cortex, and GF/SPF cecum heatmap-style figure patterns.

When a plot requires a more specialized analysis output, use this skill together with the relevant analysis skill and still preserve the publication-figure requirements here.

## Output Expectations

- Figure files and any table used for plotting.
- Plot variables, coordinate orientation, palette, and font settings.
- Any data filtering or aggregation applied only for display.
- Reused article code source, paper DOI, repository or code DOI, original file name, and dataset-specific edits.
