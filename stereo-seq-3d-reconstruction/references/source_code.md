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
