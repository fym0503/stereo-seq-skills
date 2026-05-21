#!/usr/bin/env Rscript

# Source: adapted from P09 Zenodo code file `SingleR.R`.
# Paper: Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease.
# Paper DOI: 10.1038/s41467-024-54715-y
# Code DOI: 10.5281/zenodo.14048103, license recorded by Zenodo: CC-BY-4.0.
# Reused successful pattern: SingleR reference label transfer with pruned labels,
# then export a table that can be joined back to Stereo-seq coordinates.

required_packages <- c("optparse", "Seurat", "SingleR", "scuttle", "SingleCellExperiment")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleR)
  library(scuttle)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option("--query-rds", dest = "query_rds", type = "character", help = "Seurat object for Stereo-seq bins/cellbins."),
  make_option("--reference-rds", dest = "reference_rds", type = "character", help = "Reference Seurat object."),
  make_option("--label-col", dest = "label_col", type = "character", default = "label", help = "Reference metadata column with labels."),
  make_option("--out", type = "character", default = "singler_annotation.csv", help = "Output CSV."),
  make_option("--threads", type = "integer", default = 8, help = "SingleR thread count.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$query_rds) || is.null(opt$reference_rds)) {
  stop("Required: --query-rds and --reference-rds", call. = FALSE)
}

query <- readRDS(opt$query_rds)
reference <- readRDS(opt$reference_rds)

if (!opt$label_col %in% colnames(reference@meta.data)) {
  stop(sprintf("Reference label column `%s` not found", opt$label_col), call. = FALSE)
}

reference_sce <- as.SingleCellExperiment(reference)
query_sce <- as.SingleCellExperiment(query)

reference_sce <- logNormCounts(reference_sce)
query_sce <- logNormCounts(query_sce)

prediction <- SingleR(
  test = query_sce,
  ref = reference_sce,
  labels = colData(reference_sce)[[opt$label_col]],
  de.method = "wilcox",
  num.threads = opt$threads
)

result <- data.frame(
  spatial_unit = rownames(prediction),
  label = prediction$labels,
  pruned_label = prediction$pruned.labels,
  stringsAsFactors = FALSE
)

write.csv(result, opt$out, row.names = FALSE)
message("Wrote ", opt$out)
