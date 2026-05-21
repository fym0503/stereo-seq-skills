# Source Code Registry

Use these local scripts first when CCI or ligand-receptor code is needed. Always report the reused source in the final answer.

## STcomm Mouse Organogenesis

- Paper: `Three-dimensional molecular architecture of mouse organogenesis`
- DOI: `10.1038/s41467-023-40155-7`
- Code repository: `https://github.com/gpenglab/STcomm`
- Original files: `README.md`, `R/cellColocation.R`, `R/st_comm.R`

Bundled template:

- `scripts/stcomm_colocalized_lr_template.R`: compact STcomm-style workflow for RCTD/deconvolution weights, cell-type colocalization, and CellChat LR filtering.

Reusable success points:

- Combine deconvolution-derived spatial colocalization with ligand-receptor co-expression rather than using LR expression alone.
- Keep cell-type pair colocalization table and filtered LR table as auditable outputs.
- Use Pearson/Jaccard colocalization thresholds explicitly.

## P09 AD Prefrontal Cortex CellChat

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`
- Original files: `11_1_layer_layer_communication2.ipynb`, `5_Concentric_cell_cell_communication.R`

Reusable script:

- `../stereo-seq-spatial-communication/scripts/p09_cellchat_layer_comparison.R`

## GF/SPF Cecum CellChat

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`
- Original file: `Single_cell.R`

Reusable script:

- `../stereo-seq-spatial-communication/scripts/gf_cecum_cellchat_comparison_panels.R`

## Zebrafish Heart Regeneration Spatial CellChat

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`
- Original file: `07. CellChat analysis/02.CellChat2_Stereoseq.R`

Reusable script:

- `../stereo-seq-spatial-communication/scripts/zebrafish_spatial_cellchat_template.R`

## SPIDER Spatially Variable Ligand-Receptor Interactions

- Paper: `Finding spatially variable ligand-receptor interactions with functional support from downstream genes`
- DOI: `10.1038/s41467-025-62988-0`
- Code repositories: `https://github.com/deepomicslab/SPIDER`, `https://github.com/deepomicslab/SPIDER-paper`
- Original files: `README.md`, `SPIDER-paper/notebooks/*.ipynb`

Bundled template:

- `scripts/spider_spatial_lri_template.py`: compact SPIDER wrapper for interface construction, TF scoring, spatially variable interaction testing, and plots.

Reusable success points:

- Profile ligand-receptor interfaces under spatial constraints, then test for spatially variable interaction patterns.
- Include downstream TF support when available.
- Export both interaction tables and spatial pattern/evaluation plots.

## stMLnet Multilayer CCI

- Paper: `Dissecting multilayer cell-cell communications with signaling feedback loops from spatial transcriptomics data`
- DOI: `10.1101/gr.279857.124`
- Code repository: `https://github.com/SunXQlab/stMLnet-AnalysisCode`
- Original files: `code/function.R`, `code/code.R`, `apply_in_scST/giotto_Stereoseq_dataset/s1_runstMLnet`, `benchmark/OtherMethods/NicheNet/NicheNet_main.R`, `benchmark/OtherMethods/COMMOT/commot_main.py`

Bundled template:

- `scripts/stmlnet_multilayer_cci_template.R`: parameterized multilayer CCI workflow from expression, annotation, coordinates, and prior databases.

Reusable success points:

- Model ligand-receptor, receptor-TF, and TF-target layers instead of reporting LR pairs alone.
- Calculate spatial distance matrices before pair-wise signal activity.
- Export signal activity, signal importance, and network visualizations as auditable outputs.

Environment note:

- `stMLnet-AnalysisCode` is reachable but is not a standard R package install target. If `library(stMLnet)` is unavailable, run the bundled template with `--stmlnet-source /path/to/stMLnet-AnalysisCode/code/function.R` or provide an equivalent source file.

## FlowSig Intercellular Flow

- Paper: `Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics`
- DOI: `10.1038/s41592-024-02380-w`
- Code repositories: `https://github.com/axelalmet/flowsig`, `https://github.com/axelalmet/FlowSigAnalysis_2023`
- Original files: `src/flowsig/tutorials/mouse_embryo_stereoseq_example_script.py`, `src/flowsig/preprocessing/_spatial_blocking.py`, `src/flowsig/tools/_spatial.py`, `communication_inference/analyse_spatial_communication_in_chen2022_mouse_embryo_E9.5.py`

Reusable success points:

- Use for spatially organized, pattern-driving communication flows rather than simple LR ranking.
- Keep spatial blocks, flow variables, and target gene programs explicit.
- Validate directional claims with known biology or perturbation/stage evidence.

## SpaSEG CCI

- Paper: `SpaSEG: unsupervised deep learning for multi-task analysis of spatially resolved transcriptomics`
- DOI: `10.1186/s13059-025-03697-1`
- Code repository: `https://github.com/y-bai/SpaSEG`
- Original files: `downstream/cci/cci.py`, `notebook/CCI/MB_cellChat.ipynb`, `notebook/CCI/Mouse_brain_dotplot_LR_celltype.ipynb`, `notebook/CCI/*heatmap*.ipynb`

Reusable success points:

- Use dotplots/heatmaps for LR-cell type summaries with readable labels.
- Interpret CCI together with domain and marker evidence.
- Export LR tables before figure generation.
