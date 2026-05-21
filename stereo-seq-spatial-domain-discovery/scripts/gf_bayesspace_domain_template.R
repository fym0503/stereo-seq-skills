#!/usr/bin/env Rscript

# BayesSpace spatial-domain template for Stereo-seq-like bins.
#
# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: https://github.com/1014723815/GF_SPF_cecum,
# `Spatial transcriptome.R`.
# Reusable success pattern: convert Seurat spatial counts into a
# SingleCellExperiment with coordinate colData, run BayesSpace spatialCluster,
# and export equal-aspect cluster maps with readable legends.

required_packages <- c("argparse", "Seurat", "SingleCellExperiment", "BayesSpace", "Matrix", "ggplot2", "RColorBrewer")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: BayesSpace spatial-domain discovery.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(BayesSpace)
  library(Matrix)
  library(ggplot2)
  library(RColorBrewer)
})

parser <- ArgumentParser(description = "Run BayesSpace spatial clustering on a Stereo-seq Seurat object")
parser$add_argument("--input", required = TRUE, help = "Seurat RDS")
parser$add_argument("--outdir", default = "bayesspace_domain_out")
parser$add_argument("--sample", default = "sample")
parser$add_argument("--assay", default = "Spatial")
parser$add_argument("--x-col", default = "x")
parser$add_argument("--y-col", default = "y")
parser$add_argument("--clusters", type = "integer", required = TRUE)
parser$add_argument("--pcs", type = "integer", default = 30)
parser$add_argument("--hvgs", type = "integer", default = 3000)
parser$add_argument("--nrep", type = "integer", default = 10000)
parser$add_argument("--point-size", type = "double", default = 0.4)
opts <- parser$parse_args()

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(opts$input)
counts <- as(GetAssayData(obj, assay = opts$assay, slot = "counts"), "dgCMatrix")
if (!all(c(opts$x_col, opts$y_col) %in% colnames(obj@meta.data))) {
  stop("Spatial metadata is missing coordinate columns.", call. = FALSE)
}
col_data <- data.frame(
  spot = colnames(counts),
  row = obj@meta.data[[opts$x_col]],
  col = obj@meta.data[[opts$y_col]],
  row.names = colnames(counts)
)
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene = rownames(counts), row.names = rownames(counts)),
  colData = col_data
)
set.seed(102)
sce <- spatialPreprocess(sce, platform = "ST", n.PCs = opts$pcs, n.HVGs = opts$hvgs, log.normalize = TRUE)
saveRDS(sce, file.path(opts$outdir, paste0(opts$sample, "_BayesSpace_preprocessed.rds")))
sce <- spatialCluster(
  sce,
  q = opts$clusters,
  platform = "ST",
  d = opts$pcs,
  init.method = "mclust",
  model = "t",
  gamma = 2,
  nrep = opts$nrep,
  burn.in = 200,
  save.chain = TRUE
)
saveRDS(sce, file.path(opts$outdir, paste0(opts$sample, "_BayesSpace_clustered.rds")))
df <- as.data.frame(colData(sce))
write.table(df, file.path(opts$outdir, paste0(opts$sample, "_bayesspace_clusters.tsv")), sep = "\t", quote = FALSE, row.names = TRUE)

labels <- as.character(df$spatial.cluster)
n <- length(unique(labels))
palette <- if (n <= 12) brewer.pal(max(3, n), "Set3")[seq_len(n)] else grDevices::hcl.colors(n, "Dynamic")
names(palette) <- sort(unique(labels))
p <- ggplot(df, aes(x = row, y = -col, colour = as.character(spatial.cluster))) +
  geom_point(shape = 19, size = opts$point_size, alpha = 0.95) +
  coord_fixed() +
  scale_colour_manual(values = palette) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = ifelse(n > 18, 2, 1), title = "BayesSpace")) +
  theme_void(base_size = 10, base_family = "Arial") +
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10))
ggsave(file.path(opts$outdir, paste0(opts$sample, "_BayesSpace_spatial_cluster.pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)

message("Wrote BayesSpace domain outputs to ", opts$outdir)
