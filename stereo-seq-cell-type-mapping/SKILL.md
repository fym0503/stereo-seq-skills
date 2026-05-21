---
name: stereo-seq-cell-type-mapping
description: Use when analyzing Stereo-seq or STOmics spatial transcriptomics data to annotate or map bin, spot, or cellbin units to cell types or cell states using marker genes, sc/snRNA references, RCTD, SPOTlight, cell2location, Tangram, SingleR, or related label-transfer and deconvolution workflows.
---

# Stereo-seq Cell Type Mapping

## Use This For

- Choosing a cell-type annotation, label-transfer, or deconvolution method for Stereo-seq bins, spots, or cellbins.
- Producing cell-type labels, probabilities, or abundance matrices from spatial expression data and optional sc/snRNA references.
- Auditing whether an existing Stereo-seq annotation is supported by markers, reference match, and spatial evidence.

Do not use this skill for downstream biological interpretation of mapped cell types; use `stereo-seq-spatial-cell-type-interpretation` after mapping outputs exist.

## Default Requirements

- Treat article-derived code reuse as the default when producing analysis code, plots, or tables; the user does not need to ask for this explicitly.
- Before installing packages or creating a new environment, inspect the local Python/R environments. Prefer the existing `stereo-skills-py` and `stereo-skills-r` conda environments when available: `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...`.
- If a required environment or package is missing, stop that step and tell the user exactly what is missing, which command would install it, and which analysis step is blocked. Do not silently install packages, mutate shared environments, or replace the method with an unrelated workaround unless the user asks you to.
- Before writing new plotting or analysis code, read [source_code.md](references/source_code.md) and adapt the closest bundled script in `scripts/`. Do not search GitHub/Zenodo first unless no bundled template fits.
- If no bundled article-derived script fits the task, state that explicitly and keep any custom code minimal.
- In the final response, always state which paper/code source was reused, the paper DOI, the code DOI or repository URL, the original file name, and what was changed for the current dataset.

## Workflow

1. Identify the spatial unit: DNB/bin1, bin50/bin100, pseudo-spot, segmented cellbin, or near-single-cell object.
2. Check available evidence: matched scRNA/snRNA reference, marker list, existing clusters/domains, histology, and condition/time labels.
3. Read [tool_selection.md](references/tool_selection.md) as paper-evidence context, then let the current dataset, reference availability, spatial unit, output goal, and paper-code similarity determine which template to adapt.
4. Define the input and output contract before running or reviewing code:
   - Input: expression matrix, coordinate table, spatial unit size, optional reference object, optional cluster/domain labels.
   - Output: label/probability/abundance table keyed by spatial unit, method settings, marker validation, and uncertainty flags.
5. Validate with [validation_checks.md](references/validation_checks.md).
6. If code or plotting is needed, read [source_code.md](references/source_code.md) and reuse the matching article-derived script in `scripts/` before writing new code.
7. When explaining the choice, cite the matching pilot evidence in [evidence.md](references/evidence.md), preferring Stereo-seq paper evidence over generic tool claims.

## Reusable Article Code

- `scripts/p09_singler_spatial_annotation.R`: adapted from P09 Zenodo `SingleR.R` for SingleR label transfer and `pruned.labels` export.
- `scripts/endo_seurat_label_transfer_plot.R`: adapted from Endo.R `Spatial RNA-seq analysis/01_UMAP_Label_Transfer.R` for Seurat anchor transfer and paired UMAP/spatial label plots.
- `scripts/gf_cecum_cell2location_template.py`: adapted from GF/SPF cecum `cell2location.py` for sc/snRNA signature learning and spatial abundance plots.
- `scripts/zebrafish_rctd_mapping_template.R`: adapted from zebrafish heart regeneration `04. RCTD analysis/02.RCTD.R` for spacexr/RCTD deconvolution.
- `scripts/mouse_brain_tangram_mapping_template.py`: adapted from mouse-brain Stereo-seq mapping benchmark `Code in python/tangram_circle.ipynb` for Tangram projection and probability export.

When using any bundled script, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

Always report:

- Method selected and why this route matches the data.
- Spatial unit and bin/cell size.
- Reference source and compatibility concerns.
- Main cell-type labels or abundance/probability columns.
- Markers used for validation.
- Failure modes or caveats.
- Reused article code source, including paper DOI, code DOI or repository URL, and original file name.
