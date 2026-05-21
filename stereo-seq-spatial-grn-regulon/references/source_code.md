# Source Code Registry

Use these local scripts first when GRN/regulon code or figures are needed. Always report the reused source in the final answer.

## SpaGRN

- Paper: `SpaGRN: Investigating spatially informed regulatory paths for spatially resolved transcriptomics data`
- DOI: `10.1016/j.cels.2025.101243`
- Code repository: `https://github.com/BGI-Qingdao/SpaGRN`
- Documentation/source files: `docs/source/content/01_Basic_Usage.rst`, `docs/source/Tutorials/stereo_seq_mouse_brain_hi-res.ipynb`, `docs/source/Tutorials/stereo_seq_mouse_brain_low-res.ipynb`, `src/spagrn/plot.py`, `src/spagrn/regulatory_network.py`

Bundled template:

- `scripts/spagrn_spatial_regulon_template.py`: parameterized SpaGRN inference plus top regulon spatial plots and RSS heatmaps.

Reusable success points:

- Use Stereo-seq `obsm["spatial"]` and spatial autocorrelation/coexpression during GRN inference rather than treating data as ordinary scRNA.
- Keep `n_neighbors`, autocorrelation methods, model, and raw-count layer explicit.
- Export the SpaGRN `.h5ad`, regulon dictionary, AUC matrix, and publication-style regulon plots.

## pySCENIC/SCENIC in Stereo-seq papers

- Recurrent Stereo-seq evidence: pySCENIC was observed in 25 papers in the local corpus, including mouse placentation (`10.1038/s41421-024-00740-6`) and mouse organogenesis (`10.1038/s41467-023-40155-7`).
- Representative code repository for organogenesis evidence: `https://github.com/gpenglab/STcomm`
- Tool source: `https://github.com/aertslab/pySCENIC`

Bundled template:

- `scripts/pyscenic_regulon_activity_plot_template.py`: plot pySCENIC/AUCell regulon activity matrices in spatial coordinates and cluster/domain heatmaps.
- `scripts/pyscenic_full_pipeline_template.py`: run pySCENIC GRN, motif pruning, and AUCell with explicit TF list, ranking database, and motif annotation inputs.

Reusable success points:

- Separate activity plotting from de novo GRN inference so already-computed pySCENIC outputs can be reused directly.
- Keep GRN inference, motif pruning, and AUCell outputs as separate auditable files.
- Choose top regulons by cluster/domain specificity, not only by global mean activity.
- Keep spatial maps equal-aspect and legends readable for manuscript figures.

## Wheat Root Atlas Regulons

- Paper: `A single-cell and spatial wheat root atlas with cross-species annotations delineates conserved tissue-specific marker genes and regulators`
- DOI: `10.1016/j.celrep.2025.115240`
- Code repository: `https://github.com/VIB-PSB/wheat_root_atlas`

Reusable files:

- `GRN_database/wheat MINI-ex database/01_prepare_TF2fam2mot.R`, `02_prepare_mapping.R`, `03_prepare_GOs.R`: regulon database preparation.
- `GRN_database/wheat TFs and motifs/*`: TF/motif preparation.
- `GRN_regulon_analysis/bin/cross_species_regulon_overlap.py`: cross-species regulon overlap.

Reusable success points:

- Require organism-appropriate TF, motif, and gene-motif resources.
- Separate database preparation from regulon activity analysis.
- Treat conserved regulators as hypotheses unless motif, expression, and spatial specificity agree.

## ONTraC Paper SCENIC

- Paper: `ONTraC characterizes spatially continuous variations of tissue microenvironment through niche trajectory analysis`
- DOI: `10.1186/s13059-025-03588-5`
- Code repository: `https://github.com/gyuanlab/ONTraC_paper`

Reusable files:

- `Stereo_seq_midbrain_data/post_analysis/SCENIC/Python_SCENIC.ipynb`: SCENIC analysis in a Stereo-seq midbrain workflow.
- `Stereo_seq_midbrain_data/post_analysis/Meta_cell/Meta_cell_based_gene_exp_along_trajectory.ipynb`: regulon/program interpretation along niche trajectory.

Reusable success points:

- Combine regulon activity with trajectory/domain context only after the upstream trajectory is validated.
- Export regulon activity and specificity tables before plotting.
- Keep SCENIC resource requirements explicit.
