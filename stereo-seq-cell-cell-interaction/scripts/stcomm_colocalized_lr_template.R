#!/usr/bin/env Rscript

# Template provenance:
# Paper: Three-dimensional molecular architecture of mouse organogenesis
# Paper DOI: 10.1038/s41467-023-40155-7
# Original code: gpenglab/STcomm, README.md, R/cellColocation.R, R/st_comm.R
# Reusable success pattern: combine deconvolution-derived cell-type
# colocalization with CellChat ligand-receptor evidence for spatially aware CCI.

required_packages <- c("optparse", "STcomm", "CellChat", "Seurat", "ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: STcomm colocalized LR filtering.",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(optparse)
  library(STcomm)
  library(CellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option("--seurat", type = "character", help = "Spatial Seurat object RDS."),
  make_option("--weights", type = "character", help = "CSV/TSV deconvolution/RCTD weights, spots x cell types."),
  make_option("--cellchat", type = "character", help = "CellChat object RDS or subsetCommunication table CSV/TSV."),
  make_option("--species", type = "character", default = "mouse"),
  make_option("--method", type = "character", default = "pcc"),
  make_option("--pcc", type = "double", default = 0.06),
  make_option("--jac", type = "double", default = 0.05),
  make_option("--pval", type = "double", default = 0.05),
  make_option("--padj", type = "double", default = 0.05),
  make_option("--outdir", type = "character", default = "stcomm_out")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (name in c("seurat", "weights", "cellchat")) {
  if (is.null(opt[[name]])) stop(sprintf("Required: --%s", name), call. = FALSE)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    read.delim(path, row.names = 1, check.names = FALSE)
  }
}

st_obj <- readRDS(opt$seurat)
weights <- read_table(opt$weights)
cellchat <- if (grepl("\\.rds$", opt$cellchat, ignore.case = TRUE)) readRDS(opt$cellchat) else read_table(opt$cellchat)

ctpairs <- cellColocation(
  weights,
  method = opt$method,
  pcc = opt$pcc,
  jac = opt$jac,
  pval = opt$pval,
  padj = opt$padj
)
write.csv(ctpairs, file.path(opt$outdir, "colocalized_celltype_pairs.csv"), row.names = FALSE)

st_net <- st_comm(
  object = st_obj,
  weights.df = weights,
  ctpairs = ctpairs,
  cellchat = cellchat,
  db = opt$species,
  fisher.pavl = opt$pval,
  fisher.padj = opt$padj
)
write.csv(st_net, file.path(opt$outdir, "spatially_filtered_lr_network.csv"), row.names = FALSE)

if (nrow(st_net) > 0) {
  top_net <- st_net |>
    mutate(pair = paste(source, target, sep = " -> ")) |>
    arrange(FT.padj, desc(CT.LR.pair_co.ratio)) |>
    head(30)
  p <- ggplot(top_net, aes(x = reorder(interaction_name, CT.LR.pair_co.ratio), y = CT.LR.pair_co.ratio, fill = pair)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Colocalized LR ratio", fill = "Cell-type pair") +
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  ggsave(file.path(opt$outdir, "top_spatial_lr_pairs.pdf"), p, width = 8.2, height = 6.4, dpi = 300, device = cairo_pdf)
}

message("Wrote STcomm outputs to ", opt$outdir)
