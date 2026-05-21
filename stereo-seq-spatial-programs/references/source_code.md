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
