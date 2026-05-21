#!/usr/bin/env Rscript

# stMLnet multilayer CCI template.
#
# Template provenance:
# Paper: Dissecting multilayer cell-cell communications with signaling feedback
# loops from spatial transcriptomics data
# Paper DOI: 10.1101/gr.279857.124
# Original code: https://github.com/SunXQlab/stMLnet-AnalysisCode,
# `code/function.R`, `code/code.R`, and application scripts under `apply_in_*`.
#
# Reusable success pattern: run multilayer ligand-receptor-TF-target network
# analysis from expression, annotation, coordinates, and prior databases; then
# export signal-importance tables and readable network plots.

required_packages <- c("argparse", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: stMLnet multilayer CCI.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
})

parser <- ArgumentParser(description = "Run stMLnet multilayer CCI from matrices")
parser$add_argument("--expr", required = TRUE, help = "RDS expression matrix, genes x cells")
parser$add_argument("--anno", required = TRUE, help = "RDS/data.frame annotation with columns cell and Cluster")
parser$add_argument("--coords", required = TRUE, help = "RDS/data.frame coordinates with rownames or cell column")
parser$add_argument("--database", required = TRUE, help = "RDS/RData prior database list")
parser$add_argument("--outdir", default = "stmlnet_cci_out")
parser$add_argument("--stmlnet-source", default = "", help = "Optional path to stMLnet-AnalysisCode/code/function.R or a compatible source file")
parser$add_argument("--sender", default = "", help = "Optional comma-separated sender cell types")
parser$add_argument("--receiver", default = "", help = "Optional comma-separated receiver cell types")
parser$add_argument("--ncores", type = "integer", default = 6)
parser$add_argument("--auto-para", action = "store_true", default = FALSE)
parser$add_argument("--ntrees", type = "integer", default = 500)
parser$add_argument("--ntrys", type = "integer", default = 10)
opts <- parser$parse_args()

load_stmlnet <- function(source_path) {
  if (requireNamespace("stMLnet", quietly = TRUE)) {
    suppressPackageStartupMessages(library(stMLnet))
    return(invisible(TRUE))
  }
  if (!identical(trimws(source_path), "") && file.exists(source_path)) {
    source(source_path)
    return(invisible(TRUE))
  }
  stop(
    "Missing stMLnet runtime. Install a package exposing runMLnet/getSiganlActivity/getSiganlImport ",
    "or provide --stmlnet-source pointing to stMLnet-AnalysisCode/code/function.R. ",
    "Blocked step: stMLnet multilayer CCI.",
    call. = FALSE
  )
}

load_stmlnet(opts$stmlnet_source)

read_r_object <- function(path) {
  if (grepl("\\.rda$|\\.RData$", path, ignore.case = TRUE)) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    return(env[[ls(env)[1]]])
  }
  readRDS(path)
}

split_or_null <- function(value) {
  value <- trimws(value)
  if (identical(value, "")) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
expr <- read_r_object(opts$expr)
anno <- as.data.frame(read_r_object(opts$anno))
coords <- as.data.frame(read_r_object(opts$coords))
databases <- read_r_object(opts$database)

if (!("Cluster" %in% colnames(anno))) {
  stop("Annotation table must contain a `Cluster` column.", call. = FALSE)
}
if (!("cell" %in% colnames(anno))) {
  anno$cell <- rownames(anno)
}
if (!("cell" %in% colnames(coords))) {
  coords$cell <- rownames(coords)
}
common_cells <- Reduce(intersect, list(colnames(expr), anno$cell, coords$cell))
if (!length(common_cells)) stop("No overlapping cells across expression, annotation, and coordinates.", call. = FALSE)
expr <- expr[, common_cells, drop = FALSE]
anno <- anno[match(common_cells, anno$cell), , drop = FALSE]
coords <- coords[match(common_cells, coords$cell), , drop = FALSE]
coord_numeric <- coords[, setdiff(colnames(coords), "cell"), drop = FALSE]
coord_numeric <- coord_numeric[, vapply(coord_numeric, is.numeric, logical(1)), drop = FALSE]
if (ncol(coord_numeric) < 2) stop("Coordinate table must contain at least two numeric coordinate columns.", call. = FALSE)
dist_mat <- as.matrix(dist(coord_numeric[, 1:2]))
rownames(dist_mat) <- common_cells
colnames(dist_mat) <- common_cells

sender <- split_or_null(opts$sender)
receiver <- split_or_null(opts$receiver)
res_mlnet <- runMLnet(
  ExprMat = expr,
  AnnoMat = anno,
  LigClus = sender,
  RecClus = receiver,
  Normalize = FALSE,
  OutputDir = opts$outdir,
  Databases = databases
)

mulnet_list <- list()
for (receiver_name in names(res_mlnet$mlnets)) {
  receiver_networks <- res_mlnet$mlnets[[receiver_name]]
  for (pair_name in names(receiver_networks)) {
    mlnet <- receiver_networks[[pair_name]]
    if (!is.null(mlnet$LigRec) && nrow(mlnet$LigRec) != 0) {
      mulnet_list[[pair_name]] <- mlnet
    }
  }
}
if (!length(mulnet_list)) stop("stMLnet produced no non-empty multilayer networks.", call. = FALSE)

clusters <- sort(unique(as.character(anno$Cluster)))
receiver_clusters <- if (is.null(receiver)) clusters else receiver
res_activity <- list()
for (cluster in receiver_clusters) {
  sender_clusters <- if (is.null(sender)) clusters[clusters != cluster] else sender
  res_activity[[cluster]] <- getSiganlActivity(
    ExprMat = expr,
    DistMat = dist_mat,
    AnnoMat = anno,
    MulNetList = mulnet_list,
    Receiver = cluster,
    Sender = sender_clusters,
    OutputDir = opts$outdir
  )
}
res_activity <- do.call(c, res_activity)
res_activity <- res_activity[lengths(res_activity) > 0]
if (!length(res_activity)) stop("stMLnet produced no non-empty signal activity results.", call. = FALSE)

res_importance <- list()
for (pair_name in names(res_activity)) {
  message("Calculating stMLnet signal importance: ", pair_name)
  res_importance[[pair_name]] <- getSiganlImport(
    SiganlActivity = res_activity[[pair_name]],
    Lable = pair_name,
    OutputDir = opts$outdir,
    NCores = opts$ncores,
    AutoPara = opts$auto_para,
    NTrees = opts$ntrees,
    NTrys = opts$ntrys
  )
}

result_dirs <- list(
  MLnetDir = file.path(opts$outdir, "runscMLnet"),
  ActivityDir = file.path(opts$outdir, "runModel"),
  ImportDir = file.path(opts$outdir, "getPIM")
)
saveRDS(list(mlnet = res_mlnet, activity = res_activity, importance = res_importance, dirs = result_dirs), file.path(opts$outdir, "stmlnet_results.rds"))

celltypes <- sort(unique(as.character(anno$Cluster)))
color_db <- PrepareColorDB(celltypes)
if (dir.exists(result_dirs$ImportDir)) {
  try(DrawNetworkPlot(result_dirs$ImportDir, Metric = "IM_norm", ColorDB = color_db$Celltypes, gtitle = "stMLnet CCI"), silent = TRUE)
  try(DrawHeatmapNetworkPlot(result_dirs$ImportDir, Metric = "IM_norm", ColorDB = color_db$Celltypes, gtitle = "stMLnet CCI"), silent = TRUE)
}

message("Wrote stMLnet CCI outputs to ", opts$outdir)
