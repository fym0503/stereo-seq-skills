#!/usr/bin/env Rscript

# Template provenance:
# Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the regenerating zebrafish heart
# Paper DOI: 10.1038/s41467-025-59070-0
# Original code: BGI-Qingdao/ZebrafishHeartRegeneration_project,
#   07. CellChat analysis/02.CellChat2_Stereoseq.R
# Reusable success pattern: CellChat `datatype="spatial"` with coordinates,
# distance-aware communication, per-slice circle and spatial pathway panels.

required_packages <- c("argparse", "Seurat", "CellChat")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(CellChat)
})

parser <- ArgumentParser(description = "Spatial CellChat for Stereo-seq Seurat slices")
parser$add_argument("-i", "--input", required = TRUE, help = "Seurat RDS")
parser$add_argument("-o", "--output-prefix", required = TRUE)
parser$add_argument("--split-col", default = "orig.ident")
parser$add_argument("--group-col", default = "celltype")
parser$add_argument("--x-col", default = "grid_x")
parser$add_argument("--y-col", default = "grid_y")
parser$add_argument("--assay", default = "Spatial")
parser$add_argument("--species", default = "zebrafish")
parser$add_argument("--spot-diameter", type = "double", default = 35)
parser$add_argument("--coord-scale", type = "double", default = 70)
parser$add_argument("--interaction-range", type = "double", default = 250)
parser$add_argument("--scale-distance", type = "double", default = 0.025)
parser$add_argument("--trim", type = "double", default = 0.1)
opts <- parser$parse_args()

obj <- readRDS(opts$input)
DefaultAssay(obj) <- opts$assay
obj_list <- SplitObject(obj, split.by = opts$split_col)

select_db <- function(species) {
  species <- tolower(species)
  if (species == "mouse") return(CellChatDB.mouse)
  if (species == "human") return(CellChatDB.human)
  if (exists("CellChatDB.zebrafish")) return(CellChatDB.zebrafish)
  stop("CellChatDB.zebrafish not available in this CellChat installation; provide a species-specific DB manually.", call. = FALSE)
}

for (slice in names(obj_list)) {
  seurat_object <- obj_list[[slice]]
  rdsname <- paste0(opts$output_prefix, ".", slice, ".cellchat.rds")
  if (file.exists(rdsname)) {
    cellchat <- readRDS(rdsname)
  } else {
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000, verbose = FALSE)
    data_input <- GetAssayData(seurat_object, assay = opts$assay, slot = "data")
    labels <- seurat_object@meta.data[[opts$group_col]]
    meta <- data.frame(group = labels, row.names = rownames(seurat_object@meta.data))
    spatial_loc <- data.frame(
      col = seurat_object@meta.data[[opts$x_col]] / opts$coord_scale,
      row = (seurat_object@meta.data[[opts$y_col]] * -1) / opts$coord_scale,
      row.names = rownames(seurat_object@meta.data)
    )
    scale_factors <- list(spot.diameter = opts$spot_diameter, spot = 1)
    cellchat <- createCellChat(
      object = data_input,
      meta = meta,
      group.by = "group",
      datatype = "spatial",
      coordinates = spatial_loc,
      scale.factors = scale_factors
    )
    cellchat@DB <- select_db(opts$species)
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = opts$trim,
      distance.use = TRUE,
      interaction.range = opts$interaction_range,
      scale.distance = opts$scale_distance
    )
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, file = rdsname)
  }

  group_size <- as.numeric(table(cellchat@idents))
  pdf(paste0(opts$output_prefix, ".", slice, ".circle_all.pdf"), width = 10, height = 10)
  netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
  netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  dev.off()

  pdf(paste0(opts$output_prefix, ".", slice, ".spatial_pathways.pdf"), width = 8, height = 8)
  for (path in cellchat@netP$pathways) {
    try(print(netVisual_aggregate(cellchat, layout = "spatial", signaling = path, edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)), silent = TRUE)
  }
  dev.off()
}

message("Wrote spatial CellChat panels with prefix ", opts$output_prefix)
