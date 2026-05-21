---
name: stereo-seq-quality-control-preprocessing
description: Use when Stereo-seq or STOmics data needs QC, preprocessing, GEM/GEF loading, binning, bin/cell filtering, mitochondrial/count/gene QC maps, SAW/StereoPy output handling, GEM-to-Seurat conversion, or export of cleaned h5ad/RDS objects before downstream analysis.
---

# Stereo-seq Quality Control Preprocessing

## Use This For

- Loading Stereo-seq GEM/GEF/SAW/StereoPy outputs and creating binned or filtered analysis objects.
- Plotting QC metrics in tissue coordinates before domain, mapping, CCI, trajectory, or GRN analysis.
- Converting GEM-style count coordinates into Seurat-compatible spatial objects or tabular QC summaries.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing preprocessing code or searching external repositories.
- Read [source_code.md](references/source_code.md) and match templates by input object, spatial unit, and preprocessing goal. Do not hard-code tissue-to-tool rules.
- If no curated preprocessing entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files before external search.
- Inspect local Python/R environments before running. Prefer `conda run -n stereo-skills-py python ...` for general Python scripts and `conda run -n stereo-skills-r Rscript ...` for R scripts. For StereoPy/GEF scripts that import `stereo`, use `envs/environment-python-stereopy.yml` if `stereo` is not already available. If `stereo`, `Seurat`, `Matrix`, or another required package is missing, stop and tell the user exactly what is missing and which preprocessing step is blocked.
- Preserve raw count outputs and record filtering thresholds.
- QC plots should use Arial, readable labels, equal-aspect tissue maps, and legends outside the data where possible.
- In the final response, state the reused paper, DOI, code repository or code DOI, original file, and dataset-specific edits.

## Workflow

1. Identify the raw input type: GEF, GEM, h5ad, Seurat RDS, or SAW output folder.
2. Define the spatial unit and resolution: DNB, bin, pseudo-spot, or segmented cellbin.
3. Read [source_code.md](references/source_code.md), then adapt the closest local script.
4. Export cleaned object plus QC tables and figures.
5. Report thresholds, discarded units/genes if calculated, and provenance.

## Reusable Article Code

- `scripts/endo_stereopy_qc_preprocessing_template.py`: adapted from Endo.R StereoPy binning/QC for GEF-to-h5ad preprocessing.
- `scripts/gf_gem_to_seurat_qc_template.R`: adapted from GF/SPF cecum GEM-to-Seurat and spatial QC plotting.
- `scripts/h5ad_spatial_qc_overview_template.py`: adapted from Endo.R, human cortex, and GF/SPF cecum QC/stat plotting patterns for h5ad count/feature/spatial overview figures.

## Output Expectations

- Cleaned h5ad/RDS or binned count object.
- QC tables for counts, detected genes, optional mitochondrial fraction, and coordinate bounds.
- Spatial QC maps and count/gene distributions.
- Filtering thresholds and whether they came from the template default or dataset-specific adaptation.
- Reused article code source, paper DOI, repository or code DOI, original file name, and dataset-specific edits.
