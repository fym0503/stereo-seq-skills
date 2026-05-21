# Tool Evidence Matrix

Use this as evidence for LLM reasoning, not as a fixed router. Match by biological question, object stack, desired output, and closest reusable article code.

## Spatially Variable Genes and Local Modules

Evidence pattern: the question is which genes show spatial structure or local correlation.

Evidence:
- P01 used Hotspot v0.9.1 for spatially variable genes and local-correlation gene groups.
- P10 used Hotspot v0.9.0 to identify spatially correlated modules with 2,000 HVGs and 300 neighbors.

## Spatial Co-Expression Modules

Evidence pattern: Giotto-compatible workflow and spatial co-expression module discovery objective.

Evidence:
- P08 used Giotto v3.1 `binSpect` and `clusterSpatialCorFeats` with `k=12`.

## DEG Between Groups

Evidence pattern: differential expression between domains, conditions, cell types, or time points; choose implementation by object stack and paper-code similarity.

Evidence:
- P04 and P05 used Seurat DEG.
- P09 used Scanpy `sc.tl.rank_genes_groups`.
- P10 used MAST wrapped in Seurat FindMarkers.

## Pathway or Hallmark Interpretation

Evidence pattern: gene-list enrichment or score per spatial unit, tile, or group.

Evidence:
- P04 used GSVA tumor hallmark scores.
- P10 used clusterProfiler enrichGO.

## TF/Regulon Activity

Evidence pattern: TF regulatory activity across space or domains.

Evidence:
- P02, P07, and P10 used pySCENIC-style regulon analysis and AUCell/RSS outputs.

## Latent Factors

Evidence pattern: factor-like condition-specific programs rather than individual genes or known pathways.

Evidence:
- P07 used NNMF with `rank=20`, `method="snmf/r"`, and `seed="ndsvd"` for HFD placenta section comparison.
