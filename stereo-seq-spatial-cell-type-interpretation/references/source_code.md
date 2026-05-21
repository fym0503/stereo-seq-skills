# Source Code Registry

Use these sources when adapting spatial interpretation figures. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `5_location_Tau.py`: target-cell spatial localization with highlighted target cells in red and background cells in grey.
- `case_tau_location.py` and `control_tau_location.py`: condition-stratified target-cell localization after subtype annotation.
- `2_layer_annotation.ipynb`: layer/domain proportion pies and matched spatial maps.

Reusable success points:

- Highlight focal/pathology-associated cells in red against grey background cells.
- Use the same target/background legend across raw and QC-filtered spatial maps.
- Invert coordinate axes deliberately and document the orientation.
- Use per-condition matched panels for case/control interpretation.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Single_cell.R`: condition-wise cell-type proportion summaries and comparison plots.

Bundled template:

- `scripts/gf_cecum_celltype_proportion_template.R`: stacked/fraction barplot of mapped cell types by condition or sample.

Reusable success points:

- Use stable cell-type color ordering across condition panels.
- Keep legends large enough to read in paper figures.
- Export both sample-level and condition-level proportions when sample labels are available.

## Human Cortex Spatial Cell Metadata

- Paper: `Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution`
- DOI: `10.1038/s41467-025-62793-9`
- Code repository: `https://github.com/lcy1364/Cortex-Atlas-Code`

Reusable files:

- `src/STEREO/2_Deconvolution_and_QC/4_spatialCellMeta.R`: spatial cell metadata and deconvolution interpretation plots.
- `src/STEREO/4_cellSomaProximity/somafrequentgraph_all.R`: soma/cell proximity summaries.
- `src/STEREO/4_cellSomaProximity/cellchat.R`: cell proximity and CCI context.

Reusable success points:

- Interpret cell types together with layer/domain context.
- Keep broad subclass and fine subclass interpretations separate.
- Use proximity evidence only after validating spatial mapping and coordinates.

## Avian Optic Tectum Cellbin Interpretation

- Paper: `Spatial and single-nucleus transcriptomics decoding the molecular landscape and cellular organization of avian optic tectum`
- DOI: `10.1016/j.isci.2024.109009`
- Code repository: `https://github.com/Coleliao/Spatial_OT`

Reusable files:

- `spatial_analysis/fig1BC_spatial_visualization_mannaul_annotation.R`: region/cell-type spatial annotation plots.
- `spatial_analysis/fig3BCDE_cellbin_clustering.R`: cellbin cluster interpretation.
- `spatial_analysis/fig3FGH_SGP_cellbin_subclustering.R`: subpopulation spatial interpretation.

Reusable success points:

- Interpret cell types at the same resolution used for mapping.
- Support anatomical claims with marker genes and spatial coherence.
- Report cross-species or cross-region label uncertainty explicitly.
