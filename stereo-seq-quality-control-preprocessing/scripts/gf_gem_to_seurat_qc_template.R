#!/usr/bin/env Rscript

# Stereo-seq GEM to Seurat QC template.
#
# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: https://github.com/1014723815/GF_SPF_cecum,
# `Spatial transcriptome.R`.
# Reusable success pattern: convert GEM rows to a sparse gene x coordinate
# matrix, preserve x/y coordinates in metadata, create a Seurat spatial object,
# and export readable spatial QC maps before downstream analysis.

required_packages <- c("argparse", "data.table", "Matrix", "Seurat", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: GEM to Seurat QC preprocessing.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

parser <- ArgumentParser(description = "Convert Stereo-seq GEM to Seurat object and QC maps")
parser$add_argument("--gem", required = TRUE)
parser$add_argument("--outdir", default = "gem_to_seurat_qc")
parser$add_argument("--sample", default = "sample")
parser$add_argument("--bin-size", type = "integer", default = 50)
parser$add_argument("--gene-col", default = "geneID")
parser$add_argument("--x-col", default = "x")
parser$add_argument("--y-col", default = "y")
parser$add_argument("--count-col", default = "MIDCounts")
parser$add_argument("--min-count", type = "integer", default = 100)
parser$add_argument("--min-features", type = "integer", default = 20)
opts <- parser$parse_args()

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

gem <- fread(opts$gem)
if (!all(c(opts$gene_col, opts$x_col, opts$y_col) %in% colnames(gem))) {
  stop("GEM is missing one of gene/x/y columns.", call. = FALSE)
}
if (!(opts$count_col %in% colnames(gem))) {
  fallback <- intersect(c("MIDCounts", "UMICount", "count", "counts"), colnames(gem))
  if (!length(fallback)) stop("GEM is missing count column; pass --count-col.", call. = FALSE)
  opts$count_col <- fallback[[1]]
}

setnames(gem, old = c(opts$gene_col, opts$x_col, opts$y_col, opts$count_col), new = c("geneID", "x", "y", "counts"))
gem <- gem[, .(counts = sum(counts)), by = .(geneID, x, y)]
gem[, cell := paste0(opts$sample, ":", x, "_", y)]
gene_levels <- unique(gem$geneID)
cell_levels <- unique(gem$cell)
gem[, geneIdx := match(geneID, gene_levels)]
gem[, cellIdx := match(cell, cell_levels)]

mat <- sparseMatrix(
  i = gem$geneIdx,
  j = gem$cellIdx,
  x = gem$counts,
  dimnames = list(gene_levels, cell_levels)
)
coords <- unique(gem[, .(cell, x, y)])
coords <- as.data.frame(coords)
rownames(coords) <- coords$cell

obj <- CreateSeuratObject(counts = mat, project = opts$sample, assay = "Spatial", meta.data = coords)
obj_counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
obj$nCount_Spatial <- Matrix::colSums(obj_counts)
obj$nFeature_Spatial <- Matrix::colSums(obj_counts > 0)
obj <- subset(obj, subset = nCount_Spatial >= opts$min_count & nFeature_Spatial >= opts$min_features)

saveRDS(obj, file.path(opts$outdir, paste0(opts$sample, "_bin", opts$bin_size, "_qc.rds")))
write.table(obj@meta.data, file.path(opts$outdir, paste0(opts$sample, "_qc_meta.tsv")), sep = "\t", quote = FALSE)

plot_spatial_metric <- function(meta, metric, out) {
  p <- ggplot(meta, aes(x = x, y = -y, colour = .data[[metric]])) +
    geom_point(size = 0.25, alpha = 0.9) +
    coord_fixed() +
    scale_colour_viridis_c(option = "viridis") +
    theme_void(base_size = 10, base_family = "Arial") +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(colour = metric)
  ggsave(out, plot = p, width = 6, height = 5, dpi = 300, device = cairo_pdf)
}

plot_spatial_metric(obj@meta.data, "nCount_Spatial", file.path(opts$outdir, paste0(opts$sample, "_spatial_nCount.pdf")))
plot_spatial_metric(obj@meta.data, "nFeature_Spatial", file.path(opts$outdir, paste0(opts$sample, "_spatial_nFeature.pdf")))

qc_long <- data.frame(
  metric = rep(c("nCount_Spatial", "nFeature_Spatial"), each = nrow(obj@meta.data)),
  value = c(obj$nCount_Spatial, obj$nFeature_Spatial)
)
p <- ggplot(qc_long, aes(x = metric, y = value)) +
  geom_violin(fill = "grey85", colour = "grey35", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.25) +
  theme_bw(base_size = 10, base_family = "Arial") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  ylab("Value")
ggsave(file.path(opts$outdir, paste0(opts$sample, "_qc_violin.pdf")), plot = p, width = 4, height = 4, dpi = 300, device = cairo_pdf)

message("Wrote GEM-to-Seurat QC outputs to ", opts$outdir)
