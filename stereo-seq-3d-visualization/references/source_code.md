# Source Code Registry

Use these sources before writing new 3D visualization code. Always report the reused source in the final answer.

## SPACEL Scube

- Paper: `SPACEL: deep learning-based characterization of spatial transcriptome architectures`
- DOI: `10.1038/s41467-023-43220-3`
- Code repository: `https://github.com/QuKunLab/SPACEL`
- Local corpus paper id: `S0371`

Reusable files:

- `SPACEL/Scube/plot.py`: `plot_3d`, `plot_single_slice`, and stacked slice plotting helpers.
- `SPACEL/Scube/utils_3d.py`: Open3D mesh creation, smoothing, surface sampling, color conversion, and view parameter helpers.
- `docs/tutorials/Stereo-seq_Scube.ipynb`: Stereo-seq Scube 3D alignment and visualization tutorial.
- `docs/tutorials/ST_mouse_brain_Scube.ipynb`: stacked tissue Scube workflow.

Reusable success points:

- Use Scube outputs or registered coordinates directly; do not re-register during visualization.
- Preserve slice identity and section spacing in metadata.
- Use explicit elevation/azimuth camera angles for reproducible 3D snapshots.
- For presentation panels, black backgrounds and axis-free views make sparse ST point clouds much easier to inspect.

## SLAT

- Paper: `Spatial-linked alignment tool (SLAT) for aligning heterogenous slices`
- DOI: `10.1038/s41467-023-43105-5`
- Code repository: `https://github.com/gao-lab/SLAT`
- Local corpus paper id: `S0086`

Reusable files:

- `benchmark/analysis/3d_analysis.ipynb`: 3D analysis and aligned slice visualization.
- `benchmark/analysis/plot_slices.ipynb`: slice visualization.
- `benchmark/analysis/plot_keypoints.ipynb`: keypoint visualization.

Use SLAT as upstream alignment evidence. For final rendering, pass aligned coordinates into the SPACEL/Scube visualization scripts unless reproducing SLAT figures exactly.

## STAligner

- Paper evidence: `Benchmarking clustering, alignment, and integration methods for spatial transcriptomics`
- DOI: `10.1186/s13059-024-03361-0`
- Code repository: `https://github.com/zhanglabtools/STAligner`

Reusable files:

- `Tutorials/Tutorial_3D_alignment.ipynb`: 3D alignment tutorial.

Use STAligner as upstream integration/alignment evidence. Preserve STAligner embedding/aligned-coordinate provenance when rendering the final 3D figure.

## GEM3D Toolkit

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/GEM3D_toolkit`
- Related paper code: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `gem3dtk/apply_alignment.py`: applies 3D alignment transforms to GEM/h5ad/ssDNA/cell masks.
- `gem3dtk/prepare_alignment_image.py`: prepares alignment images.
- `gemtk/draw_heatmap.py`: heatmap visualization utility.

Use GEM3D for Stereo-seq-specific coordinate generation before rendering. Keep ImageJ/TrakEM2 alignment provenance separate from visualization styling.

## Spateo

- Paper: `Spatiotemporal modeling of molecular holograms`
- Code repository: `https://github.com/aristoteleo/spateo-release`

Reusable scope:

- Use Spateo when the task requires molecular hologram, mesh, vector field, or advanced 3D embryo visualization beyond point-cloud slice snapshots.
- Treat Spateo as a heavier optional environment. Do not mix it into the SPACEL/Scube environment unless a task specifically needs it.
