#!/usr/bin/env Rscript

# Paper-quality marker/program dotplot template.
#
# Template provenance:
# Human endometrium PCOS atlas, DOI 10.1038/s41591-025-03592-z,
# https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R,
# original file `Spatial RNA-seq analysis/04.1_Marker_Visualisation_Full_Slide.R`.
# GF/SPF cecum atlas, DOI 10.1016/j.isci.2024.108941,
# https://github.com/1014723815/GF_SPF_cecum, original file
# `Spatial transcriptome.R`.
# Reusable success pattern: explicit marker/program order, readable rotated
# labels, compact legends outside the points, and vector PDF export.

required_packages <- c("argparse", "ggplot2", "readr", "dplyr", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: publication dotplot.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(scales)
})

parser <- ArgumentParser(description = "Paper-quality Stereo-seq marker/program dotplot")
parser$add_argument("--input", required = TRUE, help = "CSV/TSV long table")
parser$add_argument("--feature-col", default = "feature")
parser$add_argument("--group-col", default = "group")
parser$add_argument("--value-col", default = "avg")
parser$add_argument("--pct-col", default = "pct")
parser$add_argument("--facet-col", default = "")
parser$add_argument("--feature-order", default = "", help = "Optional comma-separated feature order")
parser$add_argument("--group-order", default = "", help = "Optional comma-separated group order")
parser$add_argument("--out", required = TRUE)
parser$add_argument("--width", type = "double", default = 8)
parser$add_argument("--height", type = "double", default = 5)
opts <- parser$parse_args()

read_any <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE)
  }
}

split_order <- function(value) {
  value <- trimws(value)
  if (identical(value, "")) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

df <- read_any(opts$input)
needed <- c(opts$feature_col, opts$group_col, opts$value_col, opts$pct_col)
missing_cols <- setdiff(needed, colnames(df))
if (length(missing_cols)) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

feature_order <- split_order(opts$feature_order)
group_order <- split_order(opts$group_order)
if (!is.null(feature_order)) df[[opts$feature_col]] <- factor(df[[opts$feature_col]], levels = feature_order)
if (!is.null(group_order)) df[[opts$group_col]] <- factor(df[[opts$group_col]], levels = group_order)

base_theme <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 10)
  )

p <- ggplot(
  df,
  aes(
    x = .data[[opts$feature_col]],
    y = .data[[opts$group_col]],
    size = .data[[opts$pct_col]],
    colour = .data[[opts$value_col]]
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(0.6, 5.2), labels = scales::label_number(accuracy = 0.1)) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
  guides(
    size = guide_legend(title = opts$pct_col, override.aes = list(colour = "grey35")),
    colour = guide_colourbar(title = opts$value_col, barheight = unit(35, "mm"))
  ) +
  base_theme

if (!identical(opts$facet_col, "") && opts$facet_col %in% colnames(df)) {
  p <- p + facet_wrap(stats::as.formula(paste("~", opts$facet_col)), scales = "free_x")
}

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
ggsave(opts$out, plot = p, width = opts$width, height = opts$height, dpi = 300, device = cairo_pdf)
message("Wrote dotplot to ", opts$out)
