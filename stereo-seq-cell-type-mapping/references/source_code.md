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

## Human Cortex Single-Cell Resolution Atlas

- Paper: `Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution`
- DOI: `10.1038/s41467-025-62793-9`
- Code repository: `https://github.com/lcy1364/Cortex-Atlas-Code`

Reusable files:

- `src/STEREO/2_Deconvolution_and_QC/1_Initial_deconvolution.R`: first-pass deconvolution.
- `src/STEREO/2_Deconvolution_and_QC/2_second_subclass_deconvolution.R`: subclass-level mapping refinement.
- `src/STEREO/2_Deconvolution_and_QC/3_RCTD_resumeMICROandOPC.R`: RCTD reruns for specific populations.
- `src/STEREO/2_Deconvolution_and_QC/4_spatialCellMeta.R`: spatial metadata and validation plots.

Reusable success points:

- Use multi-pass mapping only when the reference supports subclass resolution.
- Keep broad class, subclass, and uncertainty outputs separate.
- Validate mapped cell types with markers and spatial layer/domain context.

## Avian Optic Tectum Cross-Reference Mapping

- Paper: `Spatial and single-nucleus transcriptomics decoding the molecular landscape and cellular organization of avian optic tectum`
- DOI: `10.1016/j.isci.2024.109009`
- Code repository: `https://github.com/Coleliao/Spatial_OT`

Reusable files:

- `spatial_analysis/fig1E_integrate_spatial_sections_bin100.R`: spatial section integration.
- `spatial_analysis/fig3BCDE_cellbin_clustering.R`: cellbin clustering after segmentation.
- `snRNA_analysis/fig4D_birds_mouse_correlation.R`: cross-species/cross-reference correlation.
- `snRNA_analysis/fig4AB_single_nucleus_dOT_clustering.R`: reference clustering.

Reusable success points:

- Keep reference clustering and spatial assignment outputs auditable.
- Use correlation/orthology evidence when transferring labels across species.
- Validate spatial labels using marker genes and cellbin/bin resolution.

## SpatialID And stTransfer Evidence

- stTransfer paper: `stTransfer enables transfer of single-cell annotations to spatial transcriptomics with single-cell resolution`
- DOI: `10.1016/j.crmeth.2025.101205`
- Code repository: `https://github.com/STOmics/SpatialID`
- Benchmark evidence: `Benchmarking mapping algorithms for cell-type annotating in mouse brain by integrating single-nucleus RNA-seq and Stereo-seq data`, DOI `10.1093/bib/bbae250`

Reusable files:

- `spatialid/reader.py`, `spatialid/spatial.py`, `spatialid/trainer.py`, `spatialid/transfer.py`: SpatialID/stTransfer workflow components.
- `Code in python/spatialID.ipynb`: benchmark mapping example.

Reusable success points:

- Use when single-cell resolution or near-cellbin annotation transfer is the explicit goal.
- Export probabilities/scores rather than only hard labels.
- Compare with marker validation or a second mapping method when labels drive the main claim.

## Wheat Root Cross-Species Annotation Transfer

- Paper: `A single-cell and spatial wheat root atlas with cross-species annotations delineates conserved tissue-specific marker genes and regulators`
- DOI: `10.1016/j.celrep.2025.115240`
- Code repositories: `https://github.com/VIB-PSB/wheat_root_atlas`, `https://github.com/VIB-PSB/cross_species_annotation_transfer`

Reusable files:

- `annotation_transfer_analysis.py`: cross-species annotation transfer workflow.
- `bin/significance_deg_overlap.py`: DEG overlap significance testing.

Reusable success points:

- Require explicit ortholog mapping and gene-overlap statistics.
- Report transfer uncertainty when marker conservation is incomplete.
- Use cross-species mapping as a hypothesis generator unless independently validated.
