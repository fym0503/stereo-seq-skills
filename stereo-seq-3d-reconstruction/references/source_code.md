# Source Code Registry

Use these local scripts first when 3D reconstruction or multi-slice registration code is needed. Always report the reused source in the final answer.

## Zebrafish Heart Regeneration Atlas

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `10. Slice registrate to 3D model/CoornidateMatch_*.ipynb`: mask extraction, contour matching, affine/TPS-style coordinate transform, and registration QC.

Bundled template:

- `scripts/zebrafish_tps_slice_registration_template.py`: compact anchor/mask TPS registration and overlay plotting.

Reusable success points:

- Save transformed coordinates, not only images.
- Plot reference contours, transformed moving contours, and anchors together with equal aspect ratio.
- Prefer manually reviewed anchors when tissue shape is complex or sections are damaged.

## spatiAlign Multi-Slice Alignment

- Paper: `spatiAlign: an unsupervised contrastive learning model for data integration and alignment of spatially resolved transcriptomics`
- DOI: `10.1093/gigascience/giae042`
- Code repository: `https://github.com/STOmics/Spatialign`

Reusable files:

- `demo.py`: command-line spatiAlign training/alignment workflow.

Bundled template:

- `scripts/spatialign_multislice_alignment_template.py`: compact wrapper around the published spatiAlign workflow.

Reusable success points:

- Keep alignment parameters explicit and auditable.
- Save the trained/aligned output directory so downstream 3D reconstruction can reuse coordinates without rerunning alignment.

## GEM3D Toolkit

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/GEM3D_toolkit`

Reusable files:

- `gem3dtk/prepare_alignment_image.py`: prepare images for registration.
- `gem3dtk/apply_alignment.py`: apply registration transforms.
- `gemtk/prepare_registration_ssdna.py`: ssDNA-driven registration preparation.
- `gemRtk/STOmics_seurat.only_gem2rds.R`: GEM-to-RDS conversion support.

Reusable success points:

- Keep image preparation, registration, and coordinate transform application as separate steps.
- Save transformed GEM/cellbin coordinates for downstream analysis.
- QC registration with ssDNA/tissue-image overlays and equal-aspect plots.

## SLAT Spatial-Linked Alignment

- Paper: `Spatial-linked alignment tool (SLAT) for aligning heterogenous slices`
- DOI: `10.1038/s41467-023-43105-5`
- Code repository: `https://github.com/gao-lab/SLAT`

Reusable files:

- `benchmark/workflow/build_3d.smk`: 3D construction workflow.
- `benchmark/analysis/3d_analysis.ipynb`: 3D analysis and QC plotting.
- `benchmark/analysis/plot_slices.ipynb`, `plot_keypoints.ipynb`: slice/keypoint visualizations.
- `case/mouse_development/Seurat.ipynb`: mouse-development spatial alignment example.

Reusable success points:

- Use graph/embedding alignment only when slices/platforms need nontrivial matching.
- Export keypoints, matched pairs, and aligned coordinates.
- Validate with known anatomical landmarks before biological interpretation.

## STalign Diffeomorphic Alignment

- Paper evidence: `Benchmarking clustering, alignment, and integration methods for spatial transcriptomics`
- DOI: `10.1186/s13059-024-03361-0`
- Code repository: `https://github.com/JEFworks-Lab/STalign`

Reusable files:

- `STalign/STalign.py`: alignment implementation.
- `STalign/point_annotator.py`, `curve_annotator.py`: manual landmark/curve annotation helpers.
- `docs/notebooks/*alignment*.ipynb`: tissue-to-tissue, image-to-ST, and atlas alignment examples.

Reusable success points:

- Use landmarks/curves when automated alignment is underdetermined.
- Plot before/after overlays and deformation diagnostics.
- Keep coordinate-system orientation explicit.

## SPACEL Scube

- Paper: `SPACEL: deep learning-based characterization of spatial transcriptome architectures`
- DOI: `10.1038/s41467-023-43220-3`
- Code repository: `https://github.com/QuKunLab/SPACEL`

Reusable files:

- `docs/tutorials/Stereo-seq_Scube.ipynb`: Stereo-seq 3D architecture tutorial.
- `SPACEL/Scube/alignment.py`: slice alignment.
- `SPACEL/Scube/utils_3d.py`: 3D utilities.
- `SPACEL/Scube/plot.py`: 3D and alignment plotting helpers.

Reusable success points:

- Keep slice identities, section spacing, and transformed coordinates in metadata.
- Use 3D plots as QC plus hypothesis generation, not as proof without marker/domain validation.
- Export intermediate aligned AnnData/coordinate tables.
