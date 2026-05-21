# Source Code Registry

Use these local scripts first when trajectory code or figures are needed. Always report the reused source in the final answer.

## Zebrafish Heart Regeneration Atlas

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `06. RNA velocity analysis/03.velocity.py`: scVelo preprocessing, dynamical velocity, and velocity stream plots.
- `05. TOME analysis/04.sankey_write.R`: stage/state transition Sankey flow output.

Bundled templates:

- `scripts/zebrafish_scvelo_velocity_stream_template.py`: compact scVelo velocity workflow with metadata/UMAP import and publication-sized stream PDFs.
- `scripts/zebrafish_tome_sankey_template.R`: compact Sankey/alluvial transition plot from a transition table.

Reusable success points:

- Keep velocity plots tied to explicit embeddings and cell-state labels instead of relying on default object state.
- Use readable legend placement and Arial fonts for paper-ready output.
- Export transition links as a table alongside the Sankey HTML so the plotted flow can be audited.

## spaTrack Spatial OT Trajectory

- Paper: `Inferring cell trajectories of spatial transcriptomics via optimal transport analysis`
- DOI: `10.1016/j.cels.2025.101194`
- Code DOI: `10.5281/zenodo.14214597`
- Code repository: `https://github.com/yzf072/spaTrack`

Reusable files:

- `docs/source/notebooks/04.ST_data_of_mouse midbrain_with_multiple_times.ipynb`: multi-time ST transfer matrix and animation example.
- `spaTrack/single_time/velocity.py`: `get_ot_matrix`, `set_start_cells`, `get_ptime`, `get_velocity`.
- `spaTrack/multiple_time/transfer_matrix.py`: time-slice transition matrix utilities.

Bundled template:

- `scripts/spatrack_spatial_trajectory_template.py`: spatial OT pseudotime and trajectory stream plots from AnnData.

Reusable success points:

- Keep gene-expression and spatial-distance weights explicit.
- Require explicit start cells from cell type or coordinates.
- Export transition matrix, pseudotime table, and spatial stream plot.

## ONTraC Stereo-seq Niche Trajectory

- Paper: `ONTraC characterizes spatially continuous variations of tissue microenvironment through niche trajectory analysis`
- DOI: `10.1186/s13059-025-03588-5`
- Code DOI: `10.5281/zenodo.14171604`
- Code repositories: `https://github.com/gyuanlab/ONTraC`, `https://github.com/gyuanlab/ONTraC_paper`

Reusable files:

- `Stereo_seq_midbrain_data/run_ONTraC/stereo_midbrain_base_run_lsf.sh`: ONTraC run command for Stereo-seq mouse midbrain data.
- `Stereo_seq_midbrain_data/run_ONTraC/stereo_midbrain_base_analysis_lsf.sh`: ONTraC_analysis command and output organization.

Bundled template:

- `scripts/ontrac_stereo_niche_trajectory_template.py`: parameterized ONTraC CLI wrapper with explicit output directories and logs.

Reusable success points:

- Keep NN, GNN, and NT output directories separate.
- Preserve run logs and analysis logs.
- Treat missing ONTraC CLI as an environment blocker rather than substituting another method.
