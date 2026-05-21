#!/usr/bin/env Rscript

# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: 1014723815/GF_SPF_cecum, Single_cell.R cell-type proportion sections
# Reusable success pattern: condition-wise stacked/fraction barplot with stable
# cell-type ordering and palette, useful after cell labels are mapped spatially.

required_packages <- c("optparse", "ggplot2", "dplyr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--meta", type = "character", help = "CSV/TSV metadata table."),
  make_option("--celltype-col", dest = "celltype_col", default = "celltype"),
  make_option("--condition-col", dest = "condition_col", default = "condition"),
  make_option("--sample-col", dest = "sample_col", default = ""),
  make_option("--out", default = "celltype_proportions.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$meta)) stop("Required: --meta", call. = FALSE)

meta <- if (grepl("\\.csv$", opt$meta, ignore.case = TRUE)) read.csv(opt$meta, check.names = FALSE) else read.delim(opt$meta, check.names = FALSE)
for (col in c(opt$celltype_col, opt$condition_col)) {
  if (!col %in% colnames(meta)) stop("Missing column: ", col, call. = FALSE)
}
group_cols <- c(opt$condition_col, opt$celltype_col)
if (nzchar(opt$sample_col) && opt$sample_col %in% colnames(meta)) {
  group_cols <- c(opt$sample_col, group_cols)
}

df <- meta |>
  count(across(all_of(group_cols)), name = "n") |>
  group_by(across(all_of(setdiff(group_cols, opt$celltype_col)))) |>
  mutate(percentage = 100 * n / sum(n)) |>
  ungroup()

palette <- c(
  "#8DD3C7", "#F4F4B9", "#CFCCCF", "#D1A7B9", "#F4867C", "#86B1CD",
  "#CEB28B", "#66C2A5", "#D49A73", "#F08F6D", "#9E9DBA", "#9F9BC9",
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#57C3F3", "#476D87"
)

x_col <- if (nzchar(opt$sample_col) && opt$sample_col %in% colnames(df)) opt$sample_col else opt$condition_col
p <- ggplot(df, aes(x = .data[[x_col]], y = percentage, fill = .data[[opt$celltype_col]])) +
  geom_col(width = 0.78, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = rep(palette, length.out = length(unique(df[[opt$celltype_col]])))) +
  labs(x = NULL, y = "Cell-type proportion (%)", fill = "Cell type") +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )
if (x_col != opt$condition_col) {
  p <- p + facet_wrap(as.formula(paste("~", opt$condition_col)), ncol = 1, scales = "free_y")
}

ggsave(opt$out, p, width = 7.5, height = 5.2, dpi = 300, device = cairo_pdf)
message("Wrote ", opt$out)
