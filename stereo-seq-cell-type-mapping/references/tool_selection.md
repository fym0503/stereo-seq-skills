# Tool Evidence Matrix

Use this as evidence for LLM reasoning, not as a fixed router. Compare the user's spatial unit, reference availability, input object, output goal, and desired paper-code reuse against the entries below, then adapt the closest bundled template.

## Marker-Based Annotation

Evidence pattern: no matched scRNA/snRNA reference is available, or the task is validating cluster/domain labels.

Evidence:
- P02 annotated segmented-cell clusters with marker genes.
- P05 annotated bin clusters using markers consistent with scRNA-seq annotation.
- P08 used marker and investigated-gene evidence for tissue annotation.

Caveat: marker-only annotation is an initial or validation route; it should not be presented as high-confidence deconvolution for mixed bins.

## RCTD / spacexr

Evidence pattern: mixed bin50/bin100 spots, matched reference available, target output is cell-type composition.

Evidence:
- P01 used RCTD v1.2.0 with `doublet_mode="multi"` to resolve bin50 cell composition.
- P05 used spacexr/RCTD full mode for bin100 deconvolution from scRNA-seq reference data.

Bundled article-derived template: `scripts/zebrafish_rctd_mapping_template.R`.

## SPOTlight

Evidence pattern: R/Seurat workflow with seeded NMF-style deconvolution or paper reproduction.

Evidence:
- P03 code `single-cell-BGI/MBA/5.Spatial gene expression characterized by Stereo-seq/5.1_Stereo_unsupervised_clustering_AND_SPOTlight.R`.
- P04 used SPOTlight v0.1.6 to infer bin50 composition and assign primary cell types.

Bundled template status: evidence present in corpus; add a template only when public article code is available locally.

## cell2location

Evidence pattern: Python/anndata workflow with matched sc/snRNA reference and a goal of fine-grained abundance estimates.

Evidence:
- P07 used cell2location to map trophoblast and decidua cell types from scRNA references.
- P10 used cell2location v0.1.3 to integrate snRNA-seq and Stereo-seq.

Caveat: do not reduce the output to hard labels without preserving abundance/probability uncertainty.

## Tangram

Evidence pattern: cell/reference-to-space mapping or cross-section transfer, especially representative-section-to-section label transfer.

Evidence:
- P07 used Tangram `map_cells_to_space` to map representative subregion annotations to remaining sections.

Bundled article-derived template: `scripts/mouse_brain_tangram_mapping_template.py`.

## SingleR

Evidence pattern: cellbin or near-single-cell label transfer from reference profiles when the input unit is not heavily mixed.

Evidence:
- P09 Zenodo contains `SingleR.R` and annotation scripts.
- Original SingleR method supports reference-based annotation of single-cell profiles.

Caveat: pilot evidence is code-level and still needs manual confirmation of its exact role in P09.
