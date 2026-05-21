# Source Code Registry

Use these sources when the user asks for reusable plotting or analysis code. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Zenodo title: `Spatial Dissection of the Distinct Cellular Responses to Normal Aging and Alzheimer's Disease in Human Prefrontal Cortex at Single-Nucleus Resolution`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `SingleR.R`: reference-based layer/cell annotation using `SingleR`, `spatialLIBD`, Seurat, and `logNormCounts`.
- `3_annotation_summary.py`: Stereopy filtering, Leiden clustering, SingleR result import, and UMAP annotation plotting.
- `2_Cell_bin_ref.ipynb`: cell-type spatial scatter with fixed cell-type colors and Scanpy dotplot marker validation.

Reusable success points:

- Preserve annotation uncertainty by keeping `pruned.labels` from SingleR.
- Validate annotation with marker dotplots after label transfer.
- Use fixed categorical palettes for repeated cell types across spatial and UMAP plots.
- Keep spatial scatter equal-aspect and explicitly handle coordinate orientation.

## Human Endometrium PCOS Atlas

- Paper: `Single-cell profiling of the human endometrium in polycystic ovary syndrome`
- DOI: `10.1038/s41591-025-03592-z`
- Code repository: `https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R`

Reusable files:

- `Spatial RNA-seq analysis/01_UMAP_Label_Transfer.R`: Seurat anchor transfer from sc/snRNA reference into Stereo-seq spatial bins plus paired UMAP/spatial plots.

Bundled template:

- `scripts/endo_seurat_label_transfer_plot.R`: parameterized Seurat label-transfer workflow with PDF outputs.

Reusable success points:

- Keep reference labels explicit and export predicted labels as a table.
- Produce both reference-UMAP and spatial-coordinate views so mapping quality can be inspected.
- Use compact legends and Arial/cowplot styling for paper-sized panels.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `cell2location.py`: reference signature regression and spatial abundance mapping with cell2location.

Bundled template:

- `scripts/gf_cecum_cell2location_template.py`: parameterized cell2location reference/spatial mapping and top-abundance spatial PDFs.

Reusable success points:

- Export posterior abundance matrices for reuse by downstream interpretation/communication skills.
- Plot the dominant abundance factors directly in tissue coordinates with equal aspect ratio.
- Treat missing `cell2location` as an environment blocker instead of silently switching methods.

## Zebrafish Heart Regeneration Atlas

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `04. RCTD analysis/02.RCTD.R`: spacexr/RCTD reference construction, SpatialRNA construction from Stereo-seq coordinates, `run.RCTD(doublet_mode='full')`, and RDS export.

Bundled template:

- `scripts/zebrafish_rctd_mapping_template.R`: parameterized RCTD/spacexr mapping for Seurat spatial objects.

Reusable success points:

- Keep `Reference`, `SpatialRNA`, and RCTD object construction explicit.
- Preserve mixture weights, not only the top cell type.
- Export a top-cell-type spatial plot only as a quick QC view of the full weight matrix.

## Mouse Brain Mapping Benchmark

- Paper: `Benchmarking mapping algorithms for cell-type annotating in mouse brain by integrating single-nucleus RNA-seq and Stereo-seq data`
- DOI: `10.1093/bib/bbae250`
- Code repository: `https://github.com/qyTao185/Benchmarking-Mapping-Algorithms`
- Related method evidence: `stTransfer enables transfer of single-cell annotations to spatial transcriptomics with single-cell resolution`, DOI `10.1016/j.crmeth.2025.101205`

Reusable files:

- `Code in python/tangram_circle.ipynb`: Tangram marker selection with `rank_genes_groups`, `tg.pp_adatas`, `tg.map_cells_to_space`, `tg.project_cell_annotations`, probability export, and runtime table.

Bundled template:

- `scripts/mouse_brain_tangram_mapping_template.py`: parameterized Tangram projection and probability-map plotting.

Reusable success points:

- Select markers from the reference and explicitly intersect genes before mapping.
- Export `tangram_ct_pred` as a probability table for downstream interpretation.
- Keep runtime and marker count as audit metadata.
