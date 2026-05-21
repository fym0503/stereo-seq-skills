---
name: stereo-seq-cell-cell-interaction
description: Use when Stereo-seq or STOmics data needs cell-cell interaction, CCI, ligand-receptor, sender-receiver, niche interaction, CellChat, CellPhoneDB, NicheNet, STcomm, SPIDER, COMMOT, stMLnet, or spatially variable ligand-receptor analysis.
---

# Stereo-seq Cell-Cell Interaction

## Use This For

- CCI, ligand-receptor, sender-receiver, niche, or cell-cell communication analysis.
- CellChat/CellPhoneDB/NicheNet/STcomm/SPIDER-style workflows with spatial constraints.
- Spatially variable ligand-receptor interaction plots or spatially aware CCI filtering.

## Default Requirements

- Use bundled article-derived scripts in `scripts/` before writing new code or searching external repositories.
- Read [source_code.md](references/source_code.md) first for curated CCI templates; if no curated entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files.
- Inspect local Python/R environments first. Prefer `conda run -n stereo-skills-py python ...` for Python scripts and `conda run -n stereo-skills-r Rscript ...` for R scripts. SPIDER requires its own env (`envs/environment-python-spider.yml`) if the `spider` module is missing. If a required package is missing, stop that step and tell the user which CCI step is blocked.
- Always define upstream labels and spatial constraint before interpreting interactions.
- Keep CCI plots paper-ready: Arial, readable network/heatmap labels, no legend overlap, and vector PDF where possible.
- In the final response, state the reused paper, DOI, code repository/source file, and dataset-specific edits.

## Workflow

1. Confirm inputs: cell type/domain labels, expression object, coordinates, condition/sample labels, and optional deconvolution weights.
2. Read [source_code.md](references/source_code.md), then adapt the closest local script:
   - `scripts/stcomm_colocalized_lr_template.R` for STcomm-style colocalized cell-type and LR filtering.
   - `../stereo-seq-spatial-communication/scripts/p09_cellchat_layer_comparison.R` for P09 CellChat condition/layer comparisons.
   - `../stereo-seq-spatial-communication/scripts/gf_cecum_cellchat_comparison_panels.R` for condition-split CellChat panels.
   - `../stereo-seq-spatial-communication/scripts/zebrafish_spatial_cellchat_template.R` for spatial CellChat with coordinates.
   - `scripts/spider_spatial_lri_template.py` for SPIDER spatially variable ligand-receptor interactions.
   - `scripts/stmlnet_multilayer_cci_template.R` for stMLnet multilayer ligand-receptor-TF-target CCI with signal activity and importance.
   - `scripts/spatial_lr_neighbor_filter_template.py` for lightweight ligand-receptor support filtering by coordinate-neighbor adjacency before heavier CCI modeling.
3. Choose the template by comparing available inputs, spatial constraint, and paper-code similarity. Treat source-code/evidence files as context for LLM judgement, not as fixed tool-selection rules.
4. Validate interactions with spatial proximity, expression support, pathway context, and condition/domain specificity.
5. Report CCI results as hypotheses unless supported by spatial adjacency and orthogonal biology.

## Output Expectations

- Sender/receiver groups and spatial constraint.
- LR database/tool and species.
- Top interaction table plus network/heatmap/spatial maps.
- Condition/domain differences and caveats.
- Reused article code source, paper DOI, repository URL, original file name, and dataset-specific edits.
