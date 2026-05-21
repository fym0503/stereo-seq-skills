# Source Code Registry

Use these sources when adapting layer/domain figures. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `2_layer_annotation.ipynb`: cortical layer annotation plots, marker plots, layer proportion pies, and fixed layer colors.
- `1_bin110_analysis.ipynb`: spatial-neighbor Leiden clustering and `cluster_scatter`.
- `Graph_clustering_stereo.py`: GraphST domain clustering from Stereo-seq GEF input.
- `Layer_annotation.ipynb`: compact layer scatter entry point.

Reusable success points:

- Layer palette used in multiple P09 plots: `L1 #FF7F00`, `L2/3 #984EA3`, `L4 #4DAF4A`, `L5 #377EB8`, `L6 #E41A1C`, `WM #A65628`.
- Use ordered layer drawing so superficial/deep layer colors stay stable across spatial maps, UMAPs, and proportion plots.
- Save publication panels at `dpi=300`; hide axes for tissue maps unless coordinates are part of the claim.
- Compare domain maps with marker/QC maps such as `MBP` and `n_genes_by_counts`.

## Human Endometrium PCOS Atlas

- Paper: `Single-cell profiling of the human endometrium in polycystic ovary syndrome`
- DOI: `10.1038/s41591-025-03592-z`
- Code repository: `https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R`

Reusable files:

- `Spatial RNA-seq analysis/00.1_Stereopy_binning.py`: StereoPy GEF loading, QC, clustering, marker-gene discovery, and spatial scatter plots.

Bundled template:

- `scripts/endo_stereopy_domain_qc_template.py`: StereoPy QC/domain discovery template with marker and cluster figures.

Reusable success points:

- Keep QC, clustering, and marker plots in one auditable run directory.
- Use equal-aspect spatial scatter plots and publication DPI.
- Preserve cluster labels for later cell-type mapping and interpretation.

## Zebrafish Heart Regeneration Atlas

- Paper: `An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart`
- DOI: `10.1038/s41467-025-59070-0`
- Code repository: `https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project`

Reusable files:

- `02. Stereo-seq clustering/01.cluster.R`: multi-slice Seurat preprocessing, SCTransform, Harmony integration, clustering, and UMAP/spatial plots.

Bundled template:

- `scripts/zebrafish_harmony_domain_template.R`: parameterized multi-slice Seurat/Harmony clustering and spatial-cluster plotting.

Reusable success points:

- Keep section/sample labels in metadata before integration.
- Use Harmony only after explicit per-slice preprocessing.
- Export both integrated UMAP and per-slice spatial cluster panels.

## GraphST

- Paper: `Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST`
- DOI: `10.1038/s41467-023-36796-3`
- Code repository: `https://github.com/JinmiaoChenLab/GraphST`
- Code DOI/source record: `https://zenodo.org/record/6925603`

Reusable files:

- `GraphST/GraphST.py`: graph self-supervised representation learning with `datatype='Stereo'`.
- `GraphST/utils.py`: `clustering`, `mclust_R`, Leiden/Louvain resolution search, and optional spatial refinement.

Bundled template:

- `scripts/graphst_domain_template.py`: parameterized GraphST domain discovery for AnnData.

Reusable success points:

- Use `obsm['spatial']` explicitly and preserve the learned `obsm['emb']`.
- Export both domain labels and the trained h5ad object.
- Use equal-aspect maps with legends outside the tissue region.

## GF/SPF Cecum BayesSpace

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Spatial transcriptome.R`: BayesSpace preprocessing, spatial clustering, and cluster plotting from Stereo-seq Seurat objects.

Bundled template:

- `scripts/gf_bayesspace_domain_template.R`: parameterized BayesSpace clustering and spatial cluster PDF export.

Reusable success points:

- Keep coordinate columns and spot IDs aligned through Seurat to SingleCellExperiment conversion.
- Record `q`, `n.PCs`, `n.HVGs`, and MCMC repetition count.
- Use readable legends and equal-aspect cluster maps.

## Human Cortex Single-Cell Resolution Atlas

- Paper: `Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution`
- DOI: `10.1038/s41467-025-62793-9`
- Code repository: `https://github.com/lcy1364/Cortex-Atlas-Code`

Reusable files:

- `src/STEREO/3_domainProcess/2_generatingData.R`: domain input construction.
- `src/STEREO/3_domainProcess/DeepST_run.py`: DeepST-style spatial domain workflow.
- `src/STEREO/3_domainProcess/DeepST-main/Benchmark/run_DLPFCs_SpaGCN.py`: benchmark-style spatial domain comparison.

Reusable success points:

- Treat domain discovery and anatomical layer annotation as distinct steps.
- Use marker validation and deconvolution/QC context before naming domains.
- Keep region/layer color assignments stable across spatial maps and summaries.

## SpaSEG Domain And Marker Workflows

- Paper: `SpaSEG: unsupervised deep learning for multi-task analysis of spatially resolved transcriptomics`
- DOI: `10.1186/s13059-025-03697-1`
- Code repository: `https://github.com/y-bai/SpaSEG`

Reusable files:

- `data_processing/scanpy_processing.py`: AnnData preprocessing before model/domain analysis.
- `notebook/MultiPlatform/stereoseq_cellbin_marker_gene.ipynb`: Stereo-seq cellbin marker/domain validation.
- `downstream/plotting/_heatmapplot.py`: domain-marker heatmap plotting.

Reusable success points:

- Export domain labels and marker evidence together.
- Validate model-derived domains with known marker genes and spatial coherence.
- Use the notebook patterns as paper-derived templates while adapting all object keys locally.

## PROST Spatial Pattern And Domain Detection

- Paper: `PROST: quantitative identification of spatially variable genes and domain detection in spatial transcriptomics`
- DOI: `10.1038/s41467-024-44835-w`
- Code repository: `https://github.com/Tang-Lab-super/PROST`

Reusable files:

- `test/Stereo-seq.ipynb`: Stereo-seq example workflow.
- `PROST/plot.py`: spatial pattern/domain plotting helpers.

Reusable success points:

- Use spatially variable genes as evidence for domain labels, not only cluster boundaries.
- Keep SVG/domain outputs auditable with parameter and input-object records.
- Prefer paper-ready spatial maps with equal aspect ratio and external legends.

## Spatial Alignment And 3D-Aware Domain Methods

- SLAT paper: `Spatial-linked alignment tool (SLAT) for aligning heterogenous slices`, DOI `10.1038/s41467-023-43105-5`, repository `https://github.com/gao-lab/SLAT`
- STalign benchmark evidence: `Benchmarking clustering, alignment, and integration methods for spatial transcriptomics`, DOI `10.1186/s13059-024-03361-0`, repository `https://github.com/JEFworks-Lab/STalign`
- SPACEL paper: `SPACEL: deep learning-based characterization of spatial transcriptome architectures`, DOI `10.1038/s41467-023-43220-3`, repository `https://github.com/QuKunLab/SPACEL`

Reusable files:

- `gao-lab/SLAT`: `case/mouse_development/Seurat.ipynb`, `benchmark/analysis/plot_slices.ipynb`, `benchmark/workflow/build_3d.smk`.
- `JEFworks-Lab/STalign`: `docs/notebooks/*alignment*.ipynb`.
- `QuKunLab/SPACEL`: `docs/tutorials/Stereo-seq_Scube.ipynb`, `SPACEL/Scube/alignment.py`, `SPACEL/Scube/plot.py`.

Reusable success points:

- Use these as domain/alignment evidence only when the user has multiple slices, cross-platform sections, or atlas alignment needs.
- Keep transformed coordinates and slice/domain labels as exported tables.
- Validate alignment-derived domains against markers and tissue landmarks.
