# R Environment

Use the local conda environment for bundled R scripts:

```bash
conda run -n stereo-skills-r Rscript <script.R> ...
```

This environment is captured by `envs/environment-r-core.yml` and extended with GitHub packages in `envs/r-github-packages.R`:

- Seurat 5.3.0
- SingleR 2.4.0
- scuttle 1.12.0
- SingleCellExperiment 1.24.0
- BayesSpace 1.12.0
- spacexr/RCTD 2.2.1
- WGCNA 1.73
- hdWGCNA 0.4.11 from `smorabit/hdWGCNA`
- CellChat 2.2.0.9001 from `jinworks/CellChat`
- STcomm 0.1.2 from `gpenglab/STcomm`
- SeuratDisk from `mojaveazure/seurat-disk`
- presto 1.0.0 from `immunogenomics/presto`
- Matrix, dplyr, ggplot2, ComplexHeatmap, circlize, ggalluvial, svglite, ggpubr, networkD3, htmlwidgets

`environment-r.yml` and `envs/environment-r-core.yml` capture the conda-installable base. After recreating it, install GitHub packages:

```bash
conda run -n stereo-skills-r Rscript envs/r-github-packages.R
```

`stMLnet-AnalysisCode` is reachable but is not a standard package install target. Run `scripts/stmlnet_multilayer_cci_template.R` with `--stmlnet-source /path/to/stMLnet-AnalysisCode/code/function.R` or provide an equivalent source file.

When running a skill script, still report the article code source named in that script header.

For Python scripts, prefer `stereo-skills-py` first. It is captured by `envs/environment-python-core.yml` and includes common spatial/advanced packages such as Squidpy, decoupler, GSEApy, arboreto, scVelo, Tangram, GraphST, pySCENIC, cell2location, and spatiAlign. Use method-specific env files only for packages that intentionally conflict with the shared stack:

- SpaGRN: `envs/environment-python-spagrn.yml`
- spaTrack: `envs/environment-python-spatrack.yml`
- StereoPy: `envs/environment-python-stereopy.yml`
- SPIDER: `envs/environment-python-spider.yml`

If any required package is missing, tell the user the exact package and blocked script; do not improvise an installation unless the user asks for environment setup.
