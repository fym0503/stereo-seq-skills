---
name: stereo-seq-spatial-communication
description: Use when inferring or interpreting spatial communication in Stereo-seq data, including cell-cell interaction (CCI), ligand-receptor, bin-bin, layer-layer, interface, focal-pathology, CellChat, CellPhoneDB, NicheNet, STcomm, SPIDER, or sender-receiver analyses constrained by spatial adjacency, distance bands, tissue layers, disease regions, or hypotheses.
---

# Stereo-seq Spatial Communication

## Use This For

- Choosing communication methods for spatially resolved ligand-receptor analysis.
- Cell-cell interaction (CCI), ligand-receptor, sender-receiver, CellChat, CellPhoneDB, NicheNet, STcomm, or SPIDER workflows.
- Adding spatial constraints to CellChat, CellPhoneDB, NicheNet, or custom proximity analysis.
- Interpreting layer-layer, interface, condition, treatment, or focal-pathology communication changes.

For a task that explicitly says CCI, cell-cell interaction, ligand-receptor, STcomm, SPIDER, or spatially variable LR, also use `stereo-seq-cell-cell-interaction`.

## Default Requirements

- Treat article-derived code reuse as the default when producing analysis code, plots, or tables; the user does not need to ask for this explicitly.
- Before installing packages or creating a new environment, inspect the local Python/R environments. Prefer the existing `stereo-skills-py` and `stereo-skills-r` conda environments when available: `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...`. If using SPIDER and `spider` is missing, use `envs/environment-python-spider.yml`.
- If a required environment or package is missing, stop that step and tell the user exactly what is missing, which command would install it, and which analysis step is blocked. Do not silently install packages, mutate shared environments, or replace the method with an unrelated workaround unless the user asks you to.
- Before writing new plotting or analysis code, read [source_code.md](references/source_code.md) and adapt the closest bundled script in `scripts/`. Do not search GitHub/Zenodo first unless no bundled template fits.
- If no curated communication entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked repositories and reusable files before external search.
- If no bundled article-derived script fits the task, state that explicitly and keep any custom code minimal.
- In the final response, always state which paper/code source was reused, the paper DOI, the code DOI or repository URL, the original file name, and what was changed for the current dataset.

## Workflow

1. Confirm upstream labels: cell type, domain, layer, interface side, condition, or focal pathology.
2. Define the spatial constraint:
   - adjacent bins/cells,
   - distance threshold,
   - interface layers,
   - concentric circles,
   - same domain/layer.
3. Read [tool_selection.md](references/tool_selection.md) as paper-evidence context, then choose by task similarity, available labels, spatial constraint, and closest bundled source-code template.
4. Validate with [validation_checks.md](references/validation_checks.md).
5. If code or plotting is needed, read [source_code.md](references/source_code.md) and reuse the matching article-derived script in `scripts/` before writing new code.
6. Cite [evidence.md](references/evidence.md) when explaining the route.

## Reusable Article Code

- `scripts/p09_cellchat_layer_comparison.R`: adapted from P09 Zenodo `11_1_layer_layer_communication2.ipynb` and `5_Concentric_cell_cell_communication.R` for CellChat condition/layer comparison plots.
- `scripts/gf_cecum_cellchat_comparison_panels.R`: adapted from GF/SPF cecum `Single_cell.R` for CellChat condition-comparison circle, heatmap, and signaling-role plots.
- `scripts/zebrafish_spatial_cellchat_template.R`: adapted from zebrafish heart regeneration `07. CellChat analysis/02.CellChat2_Stereoseq.R` for spatial CellChat with coordinates and distance constraints.
- `../stereo-seq-cell-cell-interaction/scripts/stcomm_colocalized_lr_template.R`: adapted from STcomm mouse organogenesis code for deconvolution-weight colocalization plus LR filtering.
- `../stereo-seq-cell-cell-interaction/scripts/spider_spatial_lri_template.py`: adapted from SPIDER for spatially variable ligand-receptor interactions with downstream TF support.
- `scripts/spatial_lr_neighbor_filter_template.py`: adapted from STcomm, SPIDER, and P09 spatial CCI patterns for lightweight ligand-receptor support filtering by coordinate-neighbor adjacency.

When using any bundled script, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

- Sender/receiver groups.
- Ligand-receptor database and tool.
- Spatial constraint and physical units.
- Comparison design.
- Top interactions/pathways and caveats.
- Reused article code source, including paper DOI, code DOI or repository URL, and original file name.
