# Tool Evidence Matrix

Use this as evidence for LLM reasoning, not as a fixed router. Match by object type, spatial unit, domain question, and closest source-code template.

## Spatially Constrained Clustering

Evidence pattern: coordinates should influence domain discovery, especially in brain/layered tissue.

Evidence:
- P02 built a spatial KNN graph using Squidpy, combined it with expression KNN from Scanpy, then used Leiden clustering.
- P09 Zenodo includes graph clustering and layer annotation code.

## Seurat-Style Bin Clustering

Evidence pattern: R/Seurat workflow, binned matrix input, and spatial coherence validated after clustering.

Evidence:
- P03 public code creates Seurat objects from binned Stereo-seq matrices and clusters them.
- P05 and P08 Methods use Seurat normalization, scaling, clustering, and marker annotation.

## Reference-Section or Cross-Section Mapping

Evidence pattern: domains are already defined in a representative section and must be transferred to other sections.

Evidence:
- P07 used Tangram to map subregion annotations from representative sections to remaining sections.

## Manual Marker/Histology Annotation

Evidence pattern: clustering exists and the task is assigning biological domain names.

## GraphST

Evidence pattern: graph representation learning should combine gene expression and spatial coordinates before domain clustering.

Evidence:
- GraphST paper and code provide a Stereo-seq datatype path with KNN spatial graph construction, learned embeddings, clustering, and optional label refinement.

Bundled template:
- `scripts/graphst_domain_template.py`.

## BayesSpace

Evidence pattern: the workflow starts from an R/Seurat binned spatial object and seeks spatially aware cluster labels with a fixed number of domains.

Evidence:
- GF/SPF cecum public code converts Stereo-seq bins to SingleCellExperiment, runs `spatialPreprocess`, `spatialCluster`, and plots clusters in tissue coordinates.

Bundled template:
- `scripts/gf_bayesspace_domain_template.R`.

Evidence:
- P01 lung functional regions.
- P02 mouse brain anatomical regions.
- P04 liver cancer tissue regions.
- P09 cortical layers.
