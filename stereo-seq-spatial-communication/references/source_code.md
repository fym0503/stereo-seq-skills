# Source Code Registry

Use these sources when adapting communication code and figures. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `10_layer_layer_communicationR.ipynb`: layer-layer CellChat setup.
- `11_1_layer_layer_communication2.ipynb`: condition-comparison CellChat, LR strength/number plots, circle networks, heatmaps, and bubble plots.
- `5_Concentric_cell_cell_communication.R`: concentric inner/outer CellChat comparison around stress/pathology neighborhoods.

Reusable success points:

- Use `annotation` as sender/receiver identity after upstream validation.
- Use `CellChatDB.human` for human P09; switch to mouse only when species changes.
- Use `computeCommunProb(type = "triMean")`, `filterCommunication`, pathway aggregation, and `mergeCellChat` for comparison.
- Plot condition comparison with CellChat circle, heatmap, and bubble views.
- For layer plots, inherit the P09 layer palette: `L1 #FF7F00`, `L2/3 #984EA3`, `L4 #4DAF4A`, `L5 #377EB8`, `L6 #E41A1C`, `WM #A65628`.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Single_cell.R`: condition-split CellChat workflows, circle plots, heatmaps, and signaling-role comparisons.

Bundled template:

- `scripts/gf_cecum_cellchat_comparison_panels.R`: CellChat condition-comparison panels from a Seurat object and condition column.

Reusable success points:

- Split by condition only after upstream cell identities are validated.
- Export multiple complementary CellChat views instead of relying on a single network figure.
- Keep plot dimensions large enough for pathway and sender/receiver labels.

## Zebrafish Heart Regeneration Atlas

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `07. CellChat analysis/02.CellChat2_Stereoseq.R`: CellChat workflow using `datatype = "spatial"`, coordinates, distance scaling, pathway networks, and centrality plots.

Bundled template:

- `scripts/zebrafish_spatial_cellchat_template.R`: spatial CellChat template with explicit coordinate columns and distance parameters.

Reusable success points:

- Use spatial coordinates in CellChat instead of treating Stereo-seq as non-spatial scRNA.
- Keep distance scaling/spot-size parameters visible in the command line.
- Save pathway-specific network and signaling-role panels for auditability.

## STcomm Mouse Organogenesis

- Paper: `Three-dimensional molecular architecture of mouse organogenesis`
- DOI: `10.1038/s41467-023-40155-7`
- Code repository: `https://github.com/gpenglab/STcomm`
- Original files: `README.md`, `R/cellColocation.R`, `R/st_comm.R`

Reusable script:

- `../stereo-seq-cell-cell-interaction/scripts/stcomm_colocalized_lr_template.R`

## SPIDER Spatially Variable Ligand-Receptor Interactions

- Paper: `Finding spatially variable ligand-receptor interactions with functional support from downstream genes`
- DOI: `10.1038/s41467-025-62988-0`
- Code repositories: `https://github.com/deepomicslab/SPIDER`, `https://github.com/deepomicslab/SPIDER-paper`
- Original files: `README.md`, `SPIDER-paper/notebooks/*.ipynb`

Reusable script:

- `../stereo-seq-cell-cell-interaction/scripts/spider_spatial_lri_template.py`
