---
name: stereo-seq-developmental-trajectory
description: Use when Stereo-seq or STOmics data needs developmental, regeneration, temporal, pseudotime, RNA velocity, lineage-flow, cell-state transition, or stage-to-stage trajectory analysis and paper-quality trajectory plots.
---

# Stereo-seq Developmental Trajectory

## Use This For

- RNA velocity or stream plots on Stereo-seq-associated cells, bins, or matched scRNA references.
- Pseudotime, regeneration, embryonic/developmental stage, or disease-progression state transitions.
- Sankey/alluvial lineage-flow summaries between stages, domains, clusters, or inferred cell states.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing new code or searching GitHub/Zenodo.
- Read [source_code.md](references/source_code.md) first for curated trajectory templates; if no curated entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files.
- Inspect local Python/R environments first. Prefer `conda run -n stereo-skills-py python ...` for scVelo/Scanpy Python scripts and `conda run -n stereo-skills-r Rscript ...` for R scripts. spaTrack requires its own env (`envs/environment-python-spatrack.yml`) because it pins older `numpy/pandas/scipy`. If required packages are missing, stop that step and tell the user exactly which package blocks which script.
- Keep figure text paper-sized: Arial, readable legends, no legend overlap, and export PDF/300 dpi outputs where applicable.
- In the final response, state the reused paper, DOI, code repository/source file, and what was changed for the current dataset.

## Workflow

1. Identify the trajectory signal: time/stage labels, RNA velocity layers, pseudotime, lineage probability, or cluster transition table.
2. Read [source_code.md](references/source_code.md), then adapt the closest local script:
   - `scripts/zebrafish_scvelo_velocity_stream_template.py` for scVelo velocity stream plots.
   - `scripts/zebrafish_tome_sankey_template.R` for stage-to-stage Sankey flow plots.
   - `scripts/spatrack_spatial_trajectory_template.py` for spatial optimal-transport pseudotime/trajectory streams.
   - `scripts/ontrac_stereo_niche_trajectory_template.py` for ONTraC Stereo-seq niche trajectory CLI runs.
3. Let task similarity, available inputs, and paper-code provenance guide the choice; do not encode tissue-specific tool rules.
4. Preserve the input-output contract in the adapted script: cell/bin IDs, coordinates/embedding, stage labels, cell states, and output paths.
5. Validate that trajectory direction is supported by time, known markers, RNA velocity confidence, explicit start cells, or transition weights.
6. Report caveats separately from biological conclusions.

## Output Expectations

- Method and required inputs.
- Stage/state labels and trajectory direction.
- Figure files and exported tables.
- Biological interpretation with uncertainty.
- Reused article code source, paper DOI, repository URL, original file name, and dataset-specific edits.
