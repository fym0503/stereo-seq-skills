#!/usr/bin/env Rscript

# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: 1014723815/GF_SPF_cecum, Single_cell.R GO enrichment dotplot section
# Reusable success pattern: ordered enrichment dotplot with readable labels,
# horizontal coordinates, Arial base font, and compact PDF output.

required_packages <- c("optparse", "ggplot2", "dplyr", "stringr", "forcats")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(forcats)
})

option_list <- list(
  make_option("--input", type = "character", help = "CSV/TSV enrichment table."),
  make_option("--term-col", dest = "term_col", default = "Description"),
  make_option("--count-col", dest = "count_col", default = "Gene_Number"),
  make_option("--p-col", dest = "p_col", default = "p.adjust"),
  make_option("--top-n", dest = "top_n", type = "integer", default = 20),
  make_option("--out", default = "go_enrichment_dotplot.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("Required: --input", call. = FALSE)

read_table <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) read.csv(path, check.names = FALSE) else read.delim(path, check.names = FALSE)
}

df <- read_table(opt$input)
for (col in c(opt$term_col, opt$count_col, opt$p_col)) {
  if (!col %in% colnames(df)) stop("Missing column: ", col, call. = FALSE)
}
df <- df |>
  mutate(
    term = str_wrap(.data[[opt$term_col]], width = 42),
    count = as.numeric(.data[[opt$count_col]]),
    neg_log10_p = -log10(pmax(as.numeric(.data[[opt$p_col]]), .Machine$double.xmin))
  ) |>
  arrange(desc(neg_log10_p), desc(count)) |>
  slice_head(n = opt$top_n) |>
  mutate(term = fct_reorder(term, count))

height <- max(4.5, 0.28 * nrow(df) + 1.4)
p <- ggplot(df, aes(x = count, y = term, size = count, color = neg_log10_p)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#4575B4", high = "#D73027", name = "-log10(adj. p)") +
  scale_size_continuous(name = "Gene count", range = c(2.5, 7)) +
  labs(x = "Gene count", y = NULL) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = margin(8, 14, 8, 8)
  )

ggsave(opt$out, p, width = 7.2, height = height, dpi = 300, device = cairo_pdf)
message("Wrote ", opt$out)
