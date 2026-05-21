# Source Code Registry

Use these sources before writing new plotting code. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `2_layer_annotation.ipynb`: cortical layer spatial maps, marker panels, layer proportion plots, and stable ordered layer colors.
- `2_Cell_bin_ref.ipynb`: cell-type spatial scatter and marker dotplot validation.

Reusable success points:

- Ordered categorical palettes stay stable across spatial maps, UMAPs, and summary panels.
- Axes are hidden for tissue maps unless coordinates are part of the claim.
- Layer/cell-type legends are kept outside the tissue body and remain readable.

## Human Endometrium PCOS Atlas

- Paper: `Single-cell profiling of the human endometrium in polycystic ovary syndrome`
- DOI: `10.1038/s41591-025-03592-z`
- Code repository: `https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R`

Reusable files:

- `Spatial RNA-seq analysis/00.1_Stereopy_binning.py`: QC violin/spatial scatter panels from StereoPy objects.
- `Spatial RNA-seq analysis/04.1_Marker_Visualisation_Full_Slide.R`: Seurat spatial/UMAP marker dotplots, violin plots, and feature plots.

Reusable success points:

- Plot UMAP and tissue-coordinate views side by side for mapped labels.
- Use marker dotplots for annotation validation rather than only showing a spatial map.
- Keep output directories structured by sample/subset so figures remain auditable.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`
- Related cross-tissue code: `https://github.com/BGI-Intestines/Germ-free-mice`, paper DOI `10.1002/imt2.272`

Reusable files:

- `Spatial transcriptome.R`: GEM-to-Seurat conversion, BayesSpace spatial cluster maps, marker spatial plots, module-score maps.
- `Germ-free-mice/Figure*/Figure*.ipynb`: figure dotplots, violin plots, heatmaps, and condition panels.

Reusable success points:

- Equal-aspect spatial maps with cluster legends outside the tissue data.
- Dotplots use explicit group/feature order and readable axis text.
- Condition split panels use consistent palettes across GF/SPF or related contrasts.

## Human Cortex Single-Cell Resolution Atlas

- Paper: `Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution`
- DOI: `10.1038/s41467-025-62793-9`
- Code repository: `https://github.com/lcy1364/Cortex-Atlas-Code`

Reusable files:

- `src/STEREO/2_Deconvolution_and_QC/4_spatialCellMeta.R`: spatial cell metadata maps and deconvolution QC plots.
- `src/STEREO/3_domainProcess/*`: domain, layer, and marker plots after spatial clustering.
- `src/STEREO/GEM/stat.R`: GEM/bin summary plots.

Reusable success points:

- Keep cortex-region, layer, and subclass palettes stable across figures.
- Pair deconvolution/domain plots with marker and QC plots so maps are not visually unsupported.
- Use figure dimensions large enough for multi-region legends.

## Avian Optic Tectum Spatial Atlas

- Paper: `Spatial and single-nucleus transcriptomics decoding the molecular landscape and cellular organization of avian optic tectum`
- DOI: `10.1016/j.isci.2024.109009`
- Code repository: `https://github.com/Coleliao/Spatial_OT`

Reusable files:

- `functions/Visualization.R`: reusable R plotting utilities.
- `spatial_analysis/fig1BC_spatial_visualization_mannaul_annotation.R`: manual annotation spatial maps.
- `spatial_analysis/fig2D_vlnplot_selected_deg.R`: marker/DEG violin-style panels.
- `snRNA_analysis/fig4H_EXN_deg_score_on_bin100Spatial.R`: gene-score spatial maps.

Reusable success points:

- Use spatial maps, DEG/program plots, and annotation panels as a linked figure set.
- Keep layer/region labels readable and outside dense tissue areas.
- Separate bin-level and cellbin-level figures so resolution changes are explicit.

## SpaSEG Multi-Task Spatial Analysis

- Paper: `SpaSEG: unsupervised deep learning for multi-task analysis of spatially resolved transcriptomics`
- DOI: `10.1186/s13059-025-03697-1`
- Code repository: `https://github.com/y-bai/SpaSEG`

Reusable files:

- `downstream/plotting/_heatmapplot.py`: heatmap panels for markers, domains, and CCI summaries.
- `downstream/plotting/_svgplot.py`: SVG and spatial program visualization.
- `notebook/MultiPlatform/stereoseq_cellbin_marker_gene.ipynb`: Stereo-seq cellbin marker plots.

Reusable success points:

- Use one visual grammar across domain, marker, SVG, and CCI panels.
- Keep legends and heatmap labels large enough for manuscript panels.
- Use notebook outputs as figure-template evidence, but adapt paths and data objects locally.

## Thor Cell-Level Histology Plots

- Paper: `Thor: a platform for cell-level investigation of spatial transcriptomics and histology`
- DOI: `10.1038/s41467-025-62593-1`
- Code repository: `https://github.com/GuangyuWangLab2021/Thor`

Reusable files:

- `src/thor/plotting/boundary.py`: cell-boundary plotting.
- `src/thor/plotting/spot_overlap.py`: spot/bin overlap with segmented cells.
- `src/thor/plotting/spot_pie.py`: cell/spot composition pies.
- `src/thor/plotting/colors.py`: palette utilities.

Reusable success points:

- Show boundary/spot geometry when interpreting cell-level outputs.
- Avoid placing legends over histology or segmented cell boundaries.
- Use DEG and annotation plots only after multimodal alignment is inspected.
