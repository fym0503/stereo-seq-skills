# Source Code Registry

Use these sources when adapting cellbin, segmentation, mask QC, or histology-linked cell-level workflows. Always report the reused source in the final answer.

## STCellbin

- Paper: `Generating single-cell gene expression profiles for high-resolution spatial transcriptomics based on cell boundary images`
- DOI: `10.46471/gigabyte.110`
- Code repository: `https://github.com/STOmics/STCellbin`
- Original files: `STCellbin.py`, `src/cellbin/modules/cell_segmentation.py`, `src/cellbin/modules/tissue_segmentation.py`, `src/cellbin/contrib/alignment.py`, `src/cellbin/contrib/cell_mask.py`

Reusable success points:

- Treat image segmentation, tissue masking, image/expression alignment, and expression aggregation as separate auditable steps.
- Preserve model/checkpoint, mask thresholds, coordinate transform, and cell-id mapping.
- QC both mask geometry and expression counts per segmented cell before downstream cell-type mapping.

## BIDCell And CellSPA

- Paper: `BIDCell: Biologically-informed self-supervised learning for segmentation of subcellular spatial transcriptomics data`
- DOI: `10.1038/s41467-023-44560-w`
- Code repositories: `https://github.com/SydneyBioX/BIDCell`, `https://github.com/SydneyBioX/CellSPA`
- Original files: `bidcell/example_params/stereoseq.yaml`, `bidcell/processing/nuclei_segmentation.py`, `R/spatialDiversity_metrics.R`, `R/visualisation.R`

Reusable success points:

- Keep segmentation configuration in a YAML-like file so image channel, pixel scale, and biological constraints are reproducible.
- Use marker purity/diversity and negative-marker distance checks to assess whether segmentation preserves cell identity.
- Report segmentation performance metrics before making biological claims from cell-level maps.

## UCS Unified Cell Segmentation

- Paper: `UCS: A Unified Approach to Cell Segmentation for Subcellular Spatial Transcriptomics`
- DOI: `10.1002/smtd.202400975`
- Code repository: `https://github.com/YangLabHKUST/UCS`
- Related plotting/subcellular toolkit in the paper corpus: `https://github.com/ckmah/bento-tools`
- Original files from scanned code: `downstream/visualize/vizgen_mouse_brain.ipynb`, `downstream/visualize/xenium_breast_cancer.ipynb`, `bento/plotting/_plotting.py`, `bento/plotting/_layers.py`

Reusable success points:

- Use segmentation outputs as cell-level geometry, then validate with marker/spatial feature plots.
- Keep image-derived cell objects separate from expression-derived downstream clusters to avoid circular validation.
- Use compact geometry layers and readable legends for subcellular/cell-boundary figures.

## Ascidian Endostyle Stereo-seq Cell Segmentation

- Paper: `Spatially resolved single-cell atlas of ascidian endostyle provides insight into the origin of vertebrate pharyngeal organs`
- DOI: `10.1126/sciadv.adi9035`
- Code repository: `https://github.com/lskfs/ascidian-endostyle`
- Original files: `01.stereo-seq_cellseg/cellsegment.py`, `01.stereo-seq_cellseg/centroid.py`, `01.stereo-seq_cellseg/gem_mask.py`, `01.stereo-seq_cellseg/mix_segmentation.py`, `02.stereo-seq_clustering/Stereo_seurat.cellbin.R`

Reusable success points:

- Keep centroid extraction, GEM/mask intersection, mixed segmentation, and Seurat cellbin clustering as separate stages.
- Export intermediate cell coordinates and mask assignments so segmentation can be inspected independently.
- Pair cellbin clustering with spatial maps and marker validation.

## Thor Cell-Level Histology Integration

- Paper: `Thor: a platform for cell-level investigation of spatial transcriptomics and histology`
- DOI: `10.1038/s41467-025-62593-1`
- Code repositories: `https://github.com/GuangyuWangLab2021/Thor`, `https://github.com/Teichlab/bin2cell`
- Original files: `src/thor/plotting/boundary.py`, `src/thor/plotting/spot_overlap.py`, `src/thor/plotting/spot_pie.py`, `src/thor/analysis/deg.py`, `docs/_autosummary/thor.pp.nuclei_segmentation.html`, `docs/_autosummary/thor.pp.Spatial.html`

Reusable success points:

- Use boundary, spot-overlap, and spot-pie plots to show how image-derived cells relate to expression bins/spots.
- Treat histology-linked cell-level analysis as a multimodal QC problem before DEG or communication analysis.
- If Thor/bin2cell packages are unavailable, report the missing package rather than replacing cell-level analysis with ordinary bin clustering.

## Avian Optic Tectum Cellbin Analysis

- Paper: `Spatial and single-nucleus transcriptomics decoding the molecular landscape and cellular organization of avian optic tectum`
- DOI: `10.1016/j.isci.2024.109009`
- Code repository: `https://github.com/Coleliao/Spatial_OT`
- Original files: `spatial_analysis/fig3A_cell_segmentation.py`, `spatial_analysis/fig3BCDE_cellbin_clustering.R`, `spatial_analysis/fig3FGH_SGP_cellbin_subclustering.R`, `spatial_analysis/h5_to_RDS.R`

Reusable success points:

- Use cell segmentation first, then cellbin clustering/subclustering and marker validation.
- Keep section integration and cellbin-specific plots separate from bin-level spatial analysis.
- Reuse figure patterns for spatial cellbin maps with clear legends and manuscript-sized labels.
