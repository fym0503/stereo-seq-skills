---
name: stereo-seq-spatial-domain-discovery
description: Use when Stereo-seq data needs tissue domain, layer, anatomical region, disease region, or spatial cluster discovery and annotation from binned or cellbin expression matrices, coordinates, marker genes, histology, or reference-section labels.
---

# Stereo-seq Spatial Domain Discovery

## Use This For

- Discovering and naming anatomical, tissue, disease, or functional spatial domains.
- Choosing between expression-only clustering, spatially constrained clustering, Seurat-style bin clustering, or reference-section label transfer.
- Validating domains with markers, histology, spatial coherence, and differential programs.

## Default Requirements

- Treat article-derived code reuse as the default when producing analysis code, plots, or tables; the user does not need to ask for this explicitly.
- Before installing packages or creating a new environment, inspect the local Python/R environments. Prefer the existing `stereo-skills-py` and `stereo-skills-r` conda environments when available: `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...`. GraphST is configured in `stereo-skills-py`; StereoPy scripts that import `stereo` require `envs/environment-python-stereopy.yml` if `stereo` is missing.
- If a required environment or package is missing, stop that step and tell the user exactly what is missing, which command would install it, and which analysis step is blocked. Do not silently install packages, mutate shared environments, or replace the method with an unrelated workaround unless the user asks you to.
- Before writing new plotting or analysis code, read [source_code.md](references/source_code.md) and adapt the closest bundled script in `scripts/`. Do not search GitHub/Zenodo first unless no bundled template fits.
- If no curated domain entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files before external search.
- If no bundled article-derived script fits the task, state that explicitly and keep any custom code minimal.
- In the final response, always state which paper/code source was reused, the paper DOI, the code DOI or repository URL, the original file name, and what was changed for the current dataset.

## Workflow

1. Identify spatial unit and resolution: bin50/bin100/cellbin or custom pseudo-spots.
2. Decide if the task is de novo domain discovery, existing-domain annotation, or cross-section transfer.
3. Read [tool_selection.md](references/tool_selection.md) as a paper-evidence matrix, then let the current object type, spatial unit, clustering objective, and closest source-code template guide the choice.
4. Validate domain labels with [validation_checks.md](references/validation_checks.md).
5. If code or plotting is needed, read [source_code.md](references/source_code.md) and reuse the matching article-derived script in `scripts/` before writing new code.
6. Cite matching pilot evidence from [evidence.md](references/evidence.md).

## Reusable Article Code

- `scripts/p09_layer_spatial_map.py`: adapted from P09 Zenodo `2_layer_annotation.ipynb` for ordered cortical-layer spatial maps with the P09 layer palette.
- `scripts/endo_stereopy_domain_qc_template.py`: adapted from Endo.R `Spatial RNA-seq analysis/00.1_Stereopy_binning.py` for StereoPy QC, clustering, marker, and spatial-domain plots.
- `scripts/zebrafish_harmony_domain_template.R`: adapted from zebrafish heart regeneration `02. Stereo-seq clustering/01.cluster.R` for multi-slice Seurat/Harmony clustering and spatial cluster plots.
- `scripts/graphst_domain_template.py`: adapted from GraphST `GraphST/GraphST.py` and `GraphST/utils.py` for spatial graph representation learning and domain plotting.
- `scripts/gf_bayesspace_domain_template.R`: adapted from GF/SPF cecum `Spatial transcriptome.R` for BayesSpace spatial clustering and publication-ready cluster maps.
- `scripts/scanpy_spatial_domain_marker_template.py`: adapted from P09, Endo.R, and human cortex domain workflows for lightweight h5ad Leiden/domain labels, marker table, spatial map, and marker heatmap.

When using any bundled script, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

- Domain labels keyed by spatial unit.
- Method and parameters.
- Marker genes or programs supporting each domain.
- Spatial coherence and histology checks.
- Caveats for resolution, tissue mask, and cross-section transfer.
- Reused article code source, including paper DOI, code DOI or repository URL, and original file name.
