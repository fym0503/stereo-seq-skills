#!/usr/bin/env Rscript

# Paper-quality marker/program heatmap template.
#
# Template provenance:
# SpaSEG, DOI 10.1186/s13059-025-03697-1,
# https://github.com/y-bai/SpaSEG, original file
# `downstream/plotting/_heatmapplot.py` and Stereo-seq marker notebooks.
# Human cortex single-cell resolution atlas, DOI 10.1038/s41467-025-62793-9,
# https://github.com/lcy1364/Cortex-Atlas-Code, original marker/domain files
# under `src/STEREO/3_domainProcess/`.
# GF/SPF cecum atlas, DOI 10.1016/j.isci.2024.108941,
# https://github.com/1014723815/GF_SPF_cecum, original figure sections in
# `Single_cell.R`.
# Reusable success pattern: cluster/domain-by-feature matrix, readable Arial
# labels, diverging z-score color scale, and vector PDF export.

required_packages <- c("argparse", "ggplot2", "readr", "dplyr", "tidyr", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Blocked step: marker/program heatmap.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(scales)
})

parser <- ArgumentParser(description = "Paper-quality Stereo-seq marker/program heatmap")
parser$add_argument("--input", required = TRUE, help = "CSV/TSV matrix or long table.")
parser$add_argument("--format", choices = c("matrix", "long"), default = "long")
parser$add_argument("--group-col", default = "group")
parser$add_argument("--feature-col", default = "feature")
parser$add_argument("--value-col", default = "value")
parser$add_argument("--group-order", default = "")
parser$add_argument("--feature-order", default = "")
parser$add_argument("--zscore", action = "store_true", default = FALSE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--width", type = "double", default = 7)
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
if (opts$format == "matrix") {
  if (ncol(df) < 2) stop("Matrix format needs first group column plus feature columns.", call. = FALSE)
  colnames(df)[1] <- opts$group_col
  df <- df |>
    pivot_longer(cols = -all_of(opts$group_col), names_to = opts$feature_col, values_to = opts$value_col)
} else {
  needed <- c(opts$group_col, opts$feature_col, opts$value_col)
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols)) stop("Input missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

df[[opts$value_col]] <- as.numeric(df[[opts$value_col]])
if (opts$zscore) {
  df <- df |>
    group_by(.data[[opts$feature_col]]) |>
    mutate(.plot_value = as.numeric(scale(.data[[opts$value_col]]))) |>
    ungroup()
  legend_title <- "z-score"
} else {
  df$.plot_value <- df[[opts$value_col]]
  legend_title <- opts$value_col
}

group_order <- split_order(opts$group_order)
feature_order <- split_order(opts$feature_order)
if (!is.null(group_order)) df[[opts$group_col]] <- factor(df[[opts$group_col]], levels = group_order)
if (!is.null(feature_order)) df[[opts$feature_col]] <- factor(df[[opts$feature_col]], levels = feature_order)

p <- ggplot(df, aes(x = .data[[opts$feature_col]], y = .data[[opts$group_col]], fill = .plot_value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, name = legend_title) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
ggsave(opts$out, p, width = opts$width, height = opts$height, dpi = 300, device = cairo_pdf)
message("Wrote heatmap to ", opts$out)
