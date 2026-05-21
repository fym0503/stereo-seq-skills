#!/usr/bin/env Rscript

# Source: adapted from P09 Zenodo code file `9_Spatial_coexpression.R`.
# Paper: Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease.
# Paper DOI: 10.1038/s41467-024-54715-y
# Code DOI: 10.5281/zenodo.14048103, license recorded by Zenodo: CC-BY-4.0.
# Reused successful pattern: Seurat preprocessing, hdWGCNA metacells grouped by
# annotation and sample, soft-power test, module construction, hub gene output.

required_packages <- c("optparse", "Matrix", "Seurat", "WGCNA", "hdWGCNA")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
})

option_list <- list(
  make_option("--matrix", type = "character", help = "MatrixMarket file; rows are spatial units and columns are genes."),
  make_option("--meta", type = "character", help = "CSV metadata with row names or a spatial unit id column."),
  make_option("--genes", type = "character", help = "CSV with gene names in first column."),
  make_option("--group-col", dest = "group_col", type = "character", default = "annotation"),
  make_option("--sample-col", dest = "sample_col", type = "character", default = "sample"),
  make_option("--groups", type = "character", default = "", help = "Comma-separated groups to include."),
  make_option("--soft-power", dest = "soft_power", type = "integer", default = 7),
  make_option("--out-rds", dest = "out_rds", type = "character", default = "p09_hdwgcna_result.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("matrix", "meta", "genes")
for (name in required) {
  if (is.null(opt[[name]])) stop(sprintf("Required: --%s", name), call. = FALSE)
}

enableWGCNAThreads()

mat <- readMM(opt$matrix)
meta <- read.csv(opt$meta, check.names = FALSE)
genes <- read.csv(opt$genes, check.names = FALSE)[[1]]

if (!opt$group_col %in% colnames(meta)) stop("Missing group column", call. = FALSE)
if (!opt$sample_col %in% colnames(meta)) stop("Missing sample column", call. = FALSE)

rownames(mat) <- if ("X" %in% colnames(meta)) meta$X else rownames(meta)
colnames(mat) <- genes

obj <- CreateSeuratObject(counts = t(mat), min.cells = 0, min.features = 0)
for (col in colnames(meta)) obj[[col]] <- meta[[col]]

obj <- obj |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  FindNeighbors(dims = 1:30, verbose = FALSE) |>
  FindClusters(verbose = FALSE)

obj <- SetupForWGCNA(obj, gene_select = "fraction", fraction = 0.05, wgcna_name = "vis")
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c(opt$group_col, opt$sample_col),
  reduction = "pca",
  k = 25,
  max_shared = 10,
  ident.group = opt$group_col
)
obj <- NormalizeMetacells(obj)

groups <- if (nzchar(opt$groups)) strsplit(opt$groups, ",", fixed = TRUE)[[1]] else unique(obj@meta.data[[opt$group_col]])
obj <- SetDatExpr(obj, group_name = groups, group.by = opt$group_col, assay = "RNA", slot = "data")
obj <- TestSoftPowers(obj, networkType = "signed")
obj <- ConstructNetwork(obj, soft_power = opt$soft_power, setDatExpr = FALSE, tom_name = "all", overwrite_tom = TRUE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- ModuleEigengenes(obj, group.by.vars = opt$sample_col)
obj <- ModuleConnectivity(obj, group.by = opt$group_col, group_name = groups)
obj <- ResetModuleNames(obj, new_name = "Layer_modules")
obj <- ModuleExprScore(obj, n_genes = 25, method = "Seurat")

saveRDS(obj, opt$out_rds)
write.csv(GetModules(obj), sub("\\.rds$", "_modules.csv", opt$out_rds), row.names = FALSE)
write.csv(GetHubGenes(obj, n_hubs = 10), sub("\\.rds$", "_hub_genes.csv", opt$out_rds), row.names = FALSE)
message("Wrote ", opt$out_rds)
