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
