#!/usr/bin/env Rscript

# RCTD/spacexr Stereo-seq mapping template.
#
# Template provenance:
# Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart
# Paper DOI: 10.1038/s41467-025-59070-0
# Original code: https://github.com/BGI-Qingdao/ZebrafishHeartRegeneration_project,
# `04. RCTD analysis/02.RCTD.R`.
# Reusable success pattern: build a spacexr Reference from scRNA cell types,
# build SpatialRNA from Stereo-seq coordinates and counts, run RCTD, and save
# deconvolution weights for downstream spatial interpretation.

required_packages <- c("argparse", "Seurat", "spacexr", "ggplot2", "Matrix")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: RCTD/spacexr cell-type mapping.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(spacexr)
  library(ggplot2)
  library(Matrix)
})

parser <- ArgumentParser(description = "Map Stereo-seq bins/cellbins with RCTD")
parser$add_argument("--spatial", required = TRUE, help = "Spatial Seurat RDS")
parser$add_argument("--reference", required = TRUE, help = "sc/snRNA reference Seurat RDS")
parser$add_argument("--outdir", default = "rctd_mapping_out")
parser$add_argument("--sample", default = "sample")
parser$add_argument("--ref-assay", default = "RNA")
parser$add_argument("--spatial-assay", default = "Spatial")
parser$add_argument("--label-col", default = "celltype")
parser$add_argument("--exclude-label", default = "Others")
parser$add_argument("--x-col", default = "grid_x")
parser$add_argument("--y-col", default = "grid_y")
parser$add_argument("--min-count", type = "integer", default = 100)
parser$add_argument("--doublet-mode", default = "full", help = "RCTD doublet mode, e.g. full or multi")
parser$add_argument("--max-cores", type = "integer", default = 1)
opts <- parser$parse_args()

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

ref <- readRDS(opts$reference)
if (!(opts$label_col %in% colnames(ref@meta.data))) {
  stop("Reference metadata is missing label column: ", opts$label_col, call. = FALSE)
}
if (!identical(opts$exclude_label, "")) {
  ref <- subset(ref, cells = rownames(ref@meta.data)[ref@meta.data[[opts$label_col]] != opts$exclude_label])
}
ref_counts <- GetAssayData(ref, assay = opts$ref_assay, slot = "counts")
ref_celltypes <- as.factor(ref@meta.data[[opts$label_col]])
names(ref_celltypes) <- rownames(ref@meta.data)
ref_numi <- ref@meta.data[[paste0("nCount_", opts$ref_assay)]]
if (is.null(ref_numi)) ref_numi <- Matrix::colSums(ref_counts)
names(ref_numi) <- colnames(ref_counts)
reference <- Reference(ref_counts, ref_celltypes, ref_numi)

spatial <- readRDS(opts$spatial)
if (!all(c(opts$x_col, opts$y_col) %in% colnames(spatial@meta.data))) {
  stop("Spatial metadata is missing coordinate columns.", call. = FALSE)
}
count_col <- paste0("nCount_", opts$spatial_assay)
if (count_col %in% colnames(spatial@meta.data)) {
  spatial <- subset(spatial, cells = rownames(spatial@meta.data)[spatial@meta.data[[count_col]] >= opts$min_count])
}
spatial_counts <- GetAssayData(spatial, assay = opts$spatial_assay, slot = "counts")
coords <- spatial@meta.data[, c(opts$x_col, opts$y_col)]
colnames(coords) <- c("x", "y")
sp_numi <- spatial@meta.data[[count_col]]
if (is.null(sp_numi)) sp_numi <- Matrix::colSums(spatial_counts)
names(sp_numi) <- colnames(spatial_counts)
query <- SpatialRNA(coords, spatial_counts, sp_numi)

rctd <- create.RCTD(query, reference, max_cores = opts$max_cores)
rctd <- run.RCTD(rctd, doublet_mode = opts$doublet_mode)
saveRDS(rctd, file.path(opts$outdir, paste0(opts$sample, "_RCTD.rds")))

weights_mat <- as.matrix(rctd@results$weights)
weights <- as.data.frame(weights_mat)
weights$spatial_unit <- rownames(weights)
write.table(weights, file.path(opts$outdir, paste0(opts$sample, "_RCTD_weights.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

top <- weights
top$top_celltype <- colnames(weights_mat)[max.col(weights_mat, ties.method = "first")]
plot_df <- cbind(coords[match(top$spatial_unit, rownames(coords)), , drop = FALSE], top)
p <- ggplot(plot_df, aes(x = x, y = -y, colour = top_celltype)) +
  geom_point(size = 0.35, alpha = 0.9) +
  coord_fixed() +
  theme_void(base_size = 10, base_family = "Arial") +
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10)) +
  labs(colour = "Top cell type")
ggsave(file.path(opts$outdir, paste0(opts$sample, "_RCTD_top_celltype.pdf")), plot = p, width = 7, height = 5.5, dpi = 300, device = cairo_pdf)

message("Wrote RCTD mapping outputs to ", opts$outdir)
