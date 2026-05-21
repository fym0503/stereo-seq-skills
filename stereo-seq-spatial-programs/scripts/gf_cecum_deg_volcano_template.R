#!/usr/bin/env Rscript

# Template provenance:
# Paper: Effects of flora deficiency on the structure and function of the large intestine
# Paper DOI: 10.1016/j.isci.2024.108941
# Original code: 1014723815/GF_SPF_cecum, Single_cell.R DEG scatter sections
# Reusable success pattern: three-color volcano/DEG scatter with muted neutral
# points and clear up/down colors.

required_packages <- c("optparse", "ggplot2", "ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option("--input", type = "character", help = "CSV/TSV DEG table."),
  make_option("--gene-col", dest = "gene_col", default = "gene"),
  make_option("--logfc-col", dest = "logfc_col", default = "avg_log2FC"),
  make_option("--p-col", dest = "p_col", default = "p_val_adj"),
  make_option("--logfc-cutoff", dest = "logfc_cutoff", type = "double", default = 0.25),
  make_option("--p-cutoff", dest = "p_cutoff", type = "double", default = 0.05),
  make_option("--label-top", dest = "label_top", type = "integer", default = 10),
  make_option("--out", default = "deg_volcano.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("Required: --input", call. = FALSE)

df <- if (grepl("\\.csv$", opt$input, ignore.case = TRUE)) read.csv(opt$input, check.names = FALSE) else read.delim(opt$input, check.names = FALSE)
for (col in c(opt$gene_col, opt$logfc_col, opt$p_col)) {
  if (!col %in% colnames(df)) stop("Missing column: ", col, call. = FALSE)
}
df$logfc <- as.numeric(df[[opt$logfc_col]])
df$pvalue <- pmax(as.numeric(df[[opt$p_col]]), .Machine$double.xmin)
df$logP <- -log10(df$pvalue)
df$group <- "not significant"
df$group[df$logfc >= opt$logfc_cutoff & df$pvalue <= opt$p_cutoff] <- "up"
df$group[df$logfc <= -opt$logfc_cutoff & df$pvalue <= opt$p_cutoff] <- "down"

label_df <- df[df$group != "not significant", ]
label_df <- label_df[order(label_df$pvalue), ][seq_len(min(opt$label_top, nrow(label_df))), , drop = FALSE]

p <- ggplot(df, aes(x = logfc, y = logP, color = group)) +
  geom_point(size = 0.8, alpha = 0.75) +
  geom_vline(xintercept = c(-opt$logfc_cutoff, opt$logfc_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = -log10(opt$p_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = label_df, aes(label = .data[[opt$gene_col]]), size = 3, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c("down" = "#2F5688", "not significant" = "#BBBBBB", "up" = "#CC0000")) +
  labs(x = "log2 fold-change", y = "-log10 adjusted p", color = NULL) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(opt$out, p, width = 5.8, height = 4.8, dpi = 300, device = cairo_pdf)
message("Wrote ", opt$out)
