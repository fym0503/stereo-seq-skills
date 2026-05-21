#!/usr/bin/env Rscript

# Template provenance:
# Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart
# Paper DOI: 10.1038/s41467-025-59070-0
# Original code: BGI-Qingdao/ZebrafishHeartRegeneration_project,
#   02. Stereo-seq clustering/01.cluster.R
# Reusable success pattern: merge multiple Stereo-seq Seurat sections, SCTransform,
# Harmony batch correction, multi-resolution clustering, UMAP export.

required_packages <- c("argparse", "Seurat", "dplyr", "stringr", "harmony", "ggplot2", "cowplot")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(harmony)
  library(ggplot2)
  library(cowplot)
})

parser <- ArgumentParser(description = "Merge and cluster Stereo-seq Seurat slices with Harmony")
parser$add_argument("-i", "--input-dir", required = TRUE)
parser$add_argument("-o", "--outdir", required = TRUE)
parser$add_argument("-n", "--name", default = "stereo_harmony")
parser$add_argument("--assay", default = "Spatial")
parser$add_argument("--batch-col", default = "orig.ident")
parser$add_argument("--min-count", type = "integer", default = 100)
parser$add_argument("--dims", type = "integer", default = 30)
opts <- parser$parse_args()

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
rds_list <- list.files(opts$input_dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(rds_list)) stop("No .rds files found in --input-dir", call. = FALSE)

vf <- c()
objects <- list()
for (path in rds_list) {
  obj <- readRDS(path)
  DefaultAssay(obj) <- opts$assay
  obj <- subset(obj, subset = nCount_Spatial >= opts$min_count)
  obj <- SCTransform(
    obj,
    assay = opts$assay,
    variable.features.n = 2000,
    return.only.var.genes = FALSE,
    min_cells = 5,
    method = "qpoisson",
    verbose = FALSE
  )
  vf <- append(vf, VariableFeatures(obj))
  objects[[basename(path)]] <- obj
}

merged <- merge(x = objects[[1]], y = objects[2:length(objects)])
DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- unique(vf)
merged <- RunPCA(merged, verbose = FALSE)
merged <- RunHarmony(merged, group.by.vars = opts$batch_col, reduction = "pca", verbose = FALSE)
merged <- FindNeighbors(merged, dims = seq_len(opts$dims), reduction = "harmony", verbose = FALSE)
merged <- FindClusters(merged, verbose = FALSE, resolution = seq(0.1, 1.2, 0.1))
merged <- RunUMAP(merged, dims = seq_len(opts$dims), reduction = "harmony", verbose = FALSE)

saveRDS(merged, file.path(opts$outdir, paste0(opts$name, "_merge.rds")))
write.table(merged@meta.data, file.path(opts$outdir, paste0(opts$name, "_meta.tsv")), sep = "\t", quote = FALSE)

p <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "right")
ggsave(file.path(opts$outdir, paste0(opts$name, "_umap_clusters.pdf")), plot = p, width = 7, height = 5, dpi = 300, device = cairo_pdf)

if ("spatial" %in% Reductions(merged)) {
  p_sp <- DimPlot(merged, reduction = "spatial", group.by = "seurat_clusters", pt.size = 0.2) +
    theme_void(base_family = "Arial")
  ggsave(file.path(opts$outdir, paste0(opts$name, "_spatial_clusters.pdf")), plot = p_sp, width = 7, height = 6, dpi = 300, device = cairo_pdf)
}

message("Wrote Harmony domain outputs to ", opts$outdir)
