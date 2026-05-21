#!/usr/bin/env Rscript

# Template provenance:
# Paper: Single-cell profiling of the human endometrium in polycystic ovary syndrome
# Paper DOI: 10.1038/s41591-025-03592-z
# Original code: ReproductiveEndocrinologyMetabolism/Endo.R,
#   Spatial RNA-seq analysis/01_UMAP_Label_Transfer.R
# Reusable success pattern: Seurat anchor transfer from sc/snRNA reference to
# spatial bins, paired UMAP/spatial plots for predicted labels, 300 dpi PDF.

required_packages <- c("optparse", "Seurat", "SeuratDisk", "ggplot2", "cowplot")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratDisk)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--query", type = "character", help = "Spatial Seurat object (.rds or .h5seurat)."),
  make_option("--reference", type = "character", help = "Reference Seurat object (.rds or .h5seurat)."),
  make_option("--ref-label", dest = "ref_label", type = "character", default = "celltype"),
  make_option("--query-reduction", dest = "query_reduction", type = "character", default = "pca"),
  make_option("--dims", type = "character", default = "1:30"),
  make_option("--outdir", type = "character", default = "endo_label_transfer_out"),
  make_option("--prefix", type = "character", default = "spatial")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (name in c("query", "reference")) {
  if (is.null(opt[[name]])) stop(sprintf("Required: --%s", name), call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dims <- eval(parse(text = opt$dims))

load_seurat <- function(path) {
  if (grepl("\\.h5seurat$", path, ignore.case = TRUE)) {
    LoadH5Seurat(path)
  } else {
    readRDS(path)
  }
}

save_pdf <- function(plot, filename, width = 7, height = 5) {
  ggsave(filename, plot = plot, width = width, height = height, units = "in", dpi = 300, device = cairo_pdf)
}

query <- load_seurat(opt$query)
reference <- load_seurat(opt$reference)
if (!opt$ref_label %in% colnames(reference@meta.data)) {
  stop("Reference metadata does not contain --ref-label column: ", opt$ref_label, call. = FALSE)
}

DefaultAssay(reference) <- if ("RNA" %in% Assays(reference)) "RNA" else DefaultAssay(reference)
DefaultAssay(query) <- if ("RNA" %in% Assays(query)) "RNA" else DefaultAssay(query)

reference <- NormalizeData(reference, verbose = FALSE)
reference <- FindVariableFeatures(reference, verbose = FALSE)
reference <- ScaleData(reference, verbose = FALSE)
if (!"pca" %in% Reductions(reference)) {
  reference <- RunPCA(reference, verbose = FALSE)
}
if (!"umap" %in% Reductions(reference)) {
  reference <- RunUMAP(reference, dims = dims, return.model = TRUE, verbose = FALSE)
}

query <- NormalizeData(query, verbose = FALSE)
query <- FindVariableFeatures(query, verbose = FALSE)
query <- ScaleData(query, verbose = FALSE)
if (!opt$query_reduction %in% Reductions(query)) {
  query <- RunPCA(query, verbose = FALSE)
  opt$query_reduction <- "pca"
}

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = dims,
  reference.reduction = "pca"
)
ref_labels <- reference@meta.data[[opt$ref_label]]
names(ref_labels) <- colnames(reference)

query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  refdata = list(celltype = ref_labels),
  reference.reduction = "pca",
  reduction.model = "umap"
)

label_col <- "predicted.celltype"
if (!label_col %in% colnames(query@meta.data)) {
  label_col <- grep("^predicted", colnames(query@meta.data), value = TRUE)[1]
}
write.csv(query@meta.data[, label_col, drop = FALSE], file.path(opt$outdir, paste0(opt$prefix, "_predicted_labels.csv")))
saveRDS(query, file.path(opt$outdir, paste0(opt$prefix, "_label_transfer.rds")))

if ("ref.umap" %in% Reductions(query)) {
  p_umap <- DimPlot(query, reduction = "ref.umap", group.by = label_col, label = FALSE) +
    theme_cowplot(font_size = 10) +
    theme(legend.position = "right")
  save_pdf(p_umap, file.path(opt$outdir, paste0(opt$prefix, "_predicted_labels_ref_umap.pdf")))
}

spatial_reduction <- intersect(c("spatial", "Spatial", "spatial_umap"), Reductions(query))
if (length(spatial_reduction)) {
  p_spatial <- DimPlot(query, reduction = spatial_reduction[1], group.by = label_col, label = FALSE, pt.size = 0.3) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "right")
  save_pdf(p_spatial, file.path(opt$outdir, paste0(opt$prefix, "_predicted_labels_spatial.pdf")), width = 7.2, height = 6)
}

message("Wrote label-transfer results to ", opt$outdir)
