#!/usr/bin/env Rscript

# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: 1014723815/GF_SPF_cecum, Single_cell.R CellChat section
# Reusable success pattern: split Seurat cells by condition, run CellChat per
# condition, generate circle, pathway, heatmap, diff-interaction, and signaling
# role scatter panels with generous PDF dimensions.

required_packages <- c("optparse", "Seurat", "CellChat", "patchwork", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

option_list <- list(
  make_option("--seurat", type = "character", help = "Input Seurat .rds"),
  make_option("--group-col", dest = "group_col", default = "celltype"),
  make_option("--condition-col", dest = "condition_col", default = "condition"),
  make_option("--assay", default = "RNA"),
  make_option("--species", default = "mouse"),
  make_option("--pathways", default = "", help = "Comma-separated pathways to draw individually"),
  make_option("--min-cells", dest = "min_cells", type = "integer", default = 10),
  make_option("--outdir", default = "gf_cellchat_panels")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$seurat)) stop("Required: --seurat", call. = FALSE)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(opt$seurat)
DefaultAssay(obj) <- opt$assay
if (!opt$group_col %in% colnames(obj@meta.data)) stop("Missing --group-col", call. = FALSE)
if (!opt$condition_col %in% colnames(obj@meta.data)) stop("Missing --condition-col", call. = FALSE)

run_one <- function(seurat_obj, condition) {
  keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[opt$condition_col]] == condition]
  data_input <- GetAssayData(seurat_obj, assay = opt$assay, slot = "data")[, keep, drop = FALSE]
  meta <- seurat_obj@meta.data[keep, , drop = FALSE]
  meta$labels <- meta[[opt$group_col]]

  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- if (tolower(opt$species) == "mouse") CellChatDB.mouse else CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = requireNamespace("presto", quietly = TRUE))
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, if (tolower(opt$species) == "mouse") PPI.mouse else PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  saveRDS(cellchat, file.path(opt$outdir, paste0(condition, "_cellchat.rds")))
  cellchat
}

conditions <- unique(obj@meta.data[[opt$condition_col]])
objects <- setNames(lapply(conditions, function(condition) run_one(obj, condition)), conditions)

for (condition in names(objects)) {
  cellchat <- objects[[condition]]
  group_size <- as.numeric(table(cellchat@idents))
  pdf(file.path(opt$outdir, paste0(condition, "_circle_count_weight.pdf")), width = 14, height = 7)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
  dev.off()

  mat <- cellchat@net$weight
  pdf(file.path(opt$outdir, paste0(condition, "_sender_specific_circles.pdf")), width = 14, height = 10)
  n_cols <- ceiling(sqrt(nrow(mat)))
  par(mfrow = c(ceiling(nrow(mat) / n_cols), n_cols), xpd = TRUE)
  for (i in seq_len(nrow(mat))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
}

if (length(objects) >= 2) {
  merged <- mergeCellChat(objects, add.names = names(objects))
  saveRDS(merged, file.path(opt$outdir, "merged_cellchat.rds"))

  pdf(file.path(opt$outdir, "comparison_interaction_bar.pdf"), width = 6, height = 3.5)
  print(compareInteractions(merged, show.legend = FALSE, group = c(1, 2)) +
          compareInteractions(merged, show.legend = FALSE, group = c(1, 2), measure = "weight"))
  dev.off()

  pdf(file.path(opt$outdir, "comparison_diff_circle.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_diffInteraction(merged, weight.scale = TRUE)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  pdf(file.path(opt$outdir, "comparison_heatmap.pdf"), width = 9, height = 5.5)
  print(netVisual_heatmap(merged) + netVisual_heatmap(merged, measure = "weight"))
  dev.off()

  weight_minmax <- range(sapply(objects, function(x) rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)))
  role_plots <- lapply(names(objects), function(name) {
    netAnalysis_signalingRole_scatter(objects[[name]], title = name, weight.MinMax = weight_minmax)
  })
  pdf(file.path(opt$outdir, "comparison_signaling_role_scatter.pdf"), width = 11, height = 6)
  print(patchwork::wrap_plots(plots = role_plots))
  dev.off()
}

pathways <- trimws(strsplit(opt$pathways, ",", fixed = TRUE)[[1]])
pathways <- pathways[nzchar(pathways)]
for (condition in names(objects)) {
  cellchat <- objects[[condition]]
  for (pathway in pathways) {
    pdf(file.path(opt$outdir, paste0(condition, "_", pathway, "_circle.pdf")), width = 8, height = 8)
    netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    dev.off()
    pdf(file.path(opt$outdir, paste0(condition, "_", pathway, "_heatmap.pdf")), width = 7, height = 5)
    print(netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds"))
    dev.off()
  }
}

message("Wrote CellChat comparison panels to ", opt$outdir)
