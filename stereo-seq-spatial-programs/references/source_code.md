# Source Code Registry

Use these sources when adapting spatial gene-program code. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `9_Spatial_coexpression.R`: Seurat + WGCNA/hdWGCNA spatial co-expression workflow.
- `3_4_enrichR.ipynb`: enrichR pathway visualization using `plotEnrich`.
- `2_layer_annotation.ipynb`: marker dotplots and violin plots after layer annotation.

Reusable success points:

- Use Seurat normalization/PCA before hdWGCNA setup.
- Build metacells by `annotation` and `sample` to stabilize spatial co-expression modules.
- Report soft power, module names, hub genes, and module expression scores.
- Use compact enrichment plots with ordered terms and controlled label width.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Single_cell.R`: DEG volcano plots, GO/GSEA dotplots, cell-type proportion plots, and CellChat comparison panels.

Bundled templates:

- `scripts/gf_cecum_deg_volcano_template.R`: reusable DEG volcano plot with stable significance categories.
- `scripts/gf_cecum_go_gsea_dotplot_template.R`: reusable GO/GSEA dotplot with ordered terms and readable labels.

Reusable success points:

- Keep enrichment terms ordered by significance/effect size.
- Use compact, legible dotplots that avoid tiny legends and overlapping labels.
- Export vector PDF figures with Arial and 300 dpi raster fallback behavior.

## SVGbench Spatially Variable Gene Benchmark

- Paper: `Evaluating spatially variable gene detection methods for spatial transcriptomics data`
- DOI: `10.1186/s13059-023-03145-y`
- Code repository: `https://github.com/PYangLab/SVGbench`

Reusable files:

- `02_run methods/run_spatialde/mosta_spatialDE.py`: MOSTA-style SVG detection example.
- `02_run methods/run_nnSVG/mosta_nnSVG.R`: nnSVG on MOSTA-style data.
- `03_SVG detection/plot_pairwise_spearmans.R`: method concordance plots.
- `09_additional figures/plot_dataset_QC.R`: dataset QC figures for SVG benchmarking.

Reusable success points:

- Treat SVG calls as method-dependent and report thresholds/FDR settings.
- Compare top SVGs or overlap across methods when the biological claim relies on SVG robustness.
- Keep computational resource and sparsity caveats visible for large Stereo-seq data.

## PROST Spatial Pattern Genes

- Paper: `PROST: quantitative identification of spatially variable genes and domain detection in spatial transcriptomics`
- DOI: `10.1038/s41467-024-44835-w`
- Code repository: `https://github.com/Tang-Lab-super/PROST`

Reusable files:

- `test/Stereo-seq.ipynb`: Stereo-seq SVG/domain example.
- `PROST/plot.py`: spatial gene-pattern plotting helpers.

Reusable success points:

- Use SVGs to support spatial domains or anatomical gradients rather than as isolated gene lists.
- Save SVG score/rank tables before plotting top genes.
- Combine spatial expression maps with enrichment or marker context.

## Avian Optic Tectum Gene Programs

- Paper: `Spatial and single-nucleus transcriptomics decoding the molecular landscape and cellular organization of avian optic tectum`
- DOI: `10.1016/j.isci.2024.109009`
- Code repository: `https://github.com/Coleliao/Spatial_OT`

Reusable files:

- `spatial_analysis/fig2AB_differential_expression_analysis.R`: spatial DEG workflow.
- `spatial_analysis/fig2C_figS2_GO_layer_deg.R`: GO analysis for layer/domain DEGs.
- `spatial_analysis/figS3_neurotransmitter_moudle_score.R`: module-score maps.
- `snRNA_analysis/fig4EF_DEG_GO_among_dOT_EXNs.R`: DEG and GO comparison patterns.

Reusable success points:

- Link region/domain markers to GO or program summaries.
- Show spatial module scores in tissue coordinates, not only barplots.
- Keep gene-set/module names short enough for figure panels.

## Wheat Root Atlas Cross-Species Regulators

- Paper: `A single-cell and spatial wheat root atlas with cross-species annotations delineates conserved tissue-specific marker genes and regulators`
- DOI: `10.1016/j.celrep.2025.115240`
- Code repository: `https://github.com/VIB-PSB/wheat_root_atlas`

Reusable files:

- `GRN_regulon_analysis/README.md`: regulon analysis entry points.
- `GRN_regulon_analysis/bin/cross_species_regulon_overlap.py`: cross-species regulon overlap.
- `GRN_database/*`: TF, motif, OCR, and gene-motif database preparation.

Reusable success points:

- Separate marker conservation from regulator conservation.
- Require species-specific TF/motif resources before claiming conserved regulons.
- Use overlap/enrichment plots as hypothesis support, not as direct validation of target data.
