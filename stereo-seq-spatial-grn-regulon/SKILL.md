---
name: stereo-seq-spatial-grn-regulon
description: Use when Stereo-seq or STOmics spatial transcriptomics data needs gene regulatory network, transcription-factor, regulon, pySCENIC/SCENIC, SpaGRN, AUCell, regulon specificity score, or spatial TF activity analysis and paper-quality regulon plots.
---

# Stereo-seq Spatial GRN Regulon

## Use This For

- Inferring spatial gene regulatory networks from Stereo-seq expression and coordinates.
- Running or adapting SpaGRN, pySCENIC/SCENIC, AUCell, regulon specificity score, or TF activity workflows.
- Plotting regulon activity heatmaps, top regulons by cluster/domain, or regulon spatial maps.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing new code or searching external repositories.
- Inspect local Python/R environments first. Prefer `conda run -n stereo-skills-py python ...` for pySCENIC/AUCell-style Python scripts and `conda run -n stereo-skills-r Rscript ...` for R scripts. SpaGRN requires its own env (`envs/environment-python-spagrn.yml`) because it pins Python 3.8 and older `numpy/pandas/scanpy`. If a required package or database is missing, stop that step and tell the user exactly what is missing and which GRN step is blocked.
- For SpaGRN/pySCENIC, require explicit motif ranking database, motif annotation table, and TF list; do not invent those resources.
- Keep figures paper-ready: Arial, readable labels, non-overlapping legends, equal-aspect spatial plots, vector PDF when possible.
- In the final response, state the reused paper, DOI, code repository/source file, and dataset-specific edits.

## Workflow

1. Confirm input object and keys: `.h5ad`, `obsm["spatial"]`, cluster/domain label, raw-count layer, and species.
2. Read [source_code.md](references/source_code.md), then adapt the closest bundled script:
   - `scripts/spagrn_spatial_regulon_template.py` for SpaGRN spatial GRN inference and regulon activity plots.
   - `scripts/pyscenic_regulon_activity_plot_template.py` for pySCENIC/AUCell activity matrix heatmaps and spatial regulon plots.
   - `scripts/pyscenic_full_pipeline_template.py` for pySCENIC GRN, motif pruning, and AUCell pipeline execution.
3. Let the current question, available resources, and closest paper-code evidence guide the template choice; do not hard-code a tissue-specific method route.
4. Check external resources before running: TF list, `.feather` motif rankings, motif annotation `.tbl`, and optional ligand-receptor/niche table.
5. Validate outputs: number of regulons, TF-target sizes, regulon specificity scores, spatial coherence, and known marker/TF plausibility.
6. Treat GRN edges as hypotheses unless supported by motif evidence, spatial specificity, and external biology.

## Output Expectations

- Method, package versions if available, and required databases.
- Regulon table and activity matrix locations.
- Top regulons by cluster/domain and spatial maps.
- Caveats for species, motif database, sparse bins/cellbins, and small domains.
- Reused article code source, paper DOI, repository URL, original file name, and dataset-specific edits.
