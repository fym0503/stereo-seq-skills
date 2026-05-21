#!/usr/bin/env Rscript

# Source: adapted from P09 Zenodo `11_1_layer_layer_communication2.ipynb`
# and `5_Concentric_cell_cell_communication.R`.
# Paper: Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease.
# Paper DOI: 10.1038/s41467-024-54715-y
# Code DOI: 10.5281/zenodo.14048103, license recorded by Zenodo: CC-BY-4.0.
# Reused successful pattern: CellChat grouped by annotation, condition-specific
# objects, merged comparison, circle/heatmap/bubble publication panels.

required_packages <- c("optparse", "Matrix", "Seurat", "CellChat", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(CellChat)
  library(dplyr)
})

option_list <- list(
  make_option("--matrix", type = "character", help = "MatrixMarket file; rows are spatial units and columns are genes."),
  make_option("--meta", type = "character", help = "CSV metadata with annotation and condition columns."),
  make_option("--genes", type = "character", help = "CSV with gene names in first column."),
  make_option("--group-col", dest = "group_col", type = "character", default = "annotation"),
  make_option("--condition-col", dest = "condition_col", type = "character", default = "condition"),
  make_option("--species", type = "character", default = "human"),
  make_option("--min-cells", dest = "min_cells", type = "integer", default = 10),
  make_option("--outdir", type = "character", default = "p09_cellchat_out")
)

opt <- parse_args(OptionParser(option_list = option_list))
for (name in c("matrix", "meta", "genes")) {
  if (is.null(opt[[name]])) stop(sprintf("Required: --%s", name), call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

mat <- readMM(opt$matrix)
meta <- read.csv(opt$meta, check.names = FALSE)
genes <- read.csv(opt$genes, check.names = FALSE)[[1]]
rownames(mat) <- if ("X" %in% colnames(meta)) meta$X else rownames(meta)
colnames(mat) <- genes

if (!opt$group_col %in% colnames(meta)) stop("Missing group column", call. = FALSE)
if (!opt$condition_col %in% colnames(meta)) stop("Missing condition column", call. = FALSE)

data <- CreateSeuratObject(counts = t(mat), min.cells = 0, min.features = 0)
for (col in colnames(meta)) data[[col]] <- meta[[col]]
if ("samples" %in% colnames(data@meta.data)) {
  data$samples <- as.factor(data$samples)
}

run_cellchat <- function(obj) {
  count <- GetAssayData(obj, assay = "RNA", layer = "counts")
  md <- obj@meta.data
  cellchat <- createCellChat(object = count, meta = md, group.by = opt$group_col)
  cellchat <- setIdent(cellchat, ident.use = opt$group_col)
  cellchat@DB <- if (tolower(opt$species) == "mouse") CellChatDB.mouse else CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = requireNamespace("presto", quietly = TRUE))
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  aggregateNet(cellchat)
}

conditions <- unique(data@meta.data[[opt$condition_col]])
objects <- list()
for (condition in conditions) {
  keep <- rownames(data@meta.data)[data@meta.data[[opt$condition_col]] == condition]
  subset_obj <- data[, keep]
  cellchat <- run_cellchat(subset_obj)
  objects[[condition]] <- cellchat
  saveRDS(cellchat, file.path(opt$outdir, paste0(condition, "_cellchat.rds")))

  pdf(file.path(opt$outdir, paste0(condition, "_circle_count.pdf")), width = 6, height = 6)
  group_size <- as.numeric(table(cellchat@idents))
  netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = TRUE, title.name = "Number of interactions")[1]
  dev.off()
}

merged <- mergeCellChat(objects, add.names = names(objects))
saveRDS(merged, file.path(opt$outdir, "merged_cellchat.rds"))

if (length(objects) >= 2) {
  heatmap_path <- file.path(opt$outdir, "comparison_heatmap_1_vs_2.pdf")
  heatmap_error <- tryCatch(
    {
      pdf(heatmap_path, width = 4, height = 4)
      netVisual_heatmap(merged, color.heatmap = c("#2166AC", "#B2182B"), comparison = c(1, 2))
      dev.off()
      NULL
    },
    error = function(error) {
      if (dev.cur() > 1) dev.off()
      conditionMessage(error)
    }
  )
  if (!is.null(heatmap_error)) {
    if (file.exists(heatmap_path)) unlink(heatmap_path)
    writeLines(heatmap_error, file.path(opt$outdir, "comparison_heatmap_1_vs_2.error.txt"))
    warning("CellChat comparison heatmap failed; CellChat objects and circle plots were still written.")
  }
}

message("Wrote CellChat objects and plots to ", opt$outdir)
