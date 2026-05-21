#!/usr/bin/env Rscript

# Template provenance:
# Paper: An organ-wide spatiotemporal transcriptomic and cellular atlas of the
# regenerating zebrafish heart
# Paper DOI: 10.1038/s41467-025-59070-0
# Original code: BGI-Qingdao/ZebrafishHeartRegeneration_project,
#   05. TOME analysis/04.sankey_write.R
# Reusable success pattern: convert stage/state transitions into auditable
# source-target-value links and render a Sankey figure.

required_packages <- c("optparse", "dplyr", "readr", "networkD3", "htmlwidgets")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(networkD3)
  library(htmlwidgets)
})

option_list <- list(
  make_option("--links", type = "character", help = "CSV/TSV with source,target,value columns."),
  make_option("--source-col", dest = "source_col", default = "source"),
  make_option("--target-col", dest = "target_col", default = "target"),
  make_option("--value-col", dest = "value_col", default = "value"),
  make_option("--min-value", dest = "min_value", type = "double", default = 0),
  make_option("--out-html", dest = "out_html", default = "trajectory_sankey.html"),
  make_option("--out-links", dest = "out_links", default = "trajectory_sankey_links.tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$links)) stop("Required: --links", call. = FALSE)

links_raw <- read_delim(opt$links, delim = ifelse(grepl("\\.csv$", opt$links, ignore.case = TRUE), ",", "\t"), show_col_types = FALSE)
for (col in c(opt$source_col, opt$target_col, opt$value_col)) {
  if (!col %in% colnames(links_raw)) stop("Missing column: ", col, call. = FALSE)
}

links <- links_raw |>
  transmute(
    source = as.character(.data[[opt$source_col]]),
    target = as.character(.data[[opt$target_col]]),
    value = as.numeric(.data[[opt$value_col]])
  ) |>
  filter(!is.na(source), !is.na(target), !is.na(value), value > opt$min_value) |>
  group_by(source, target) |>
  summarise(value = sum(value), .groups = "drop")

nodes <- data.frame(name = unique(c(links$source, links$target)), stringsAsFactors = FALSE)
links_plot <- links |>
  mutate(
    source_id = match(source, nodes$name) - 1,
    target_id = match(target, nodes$name) - 1
  )

write_tsv(links, opt$out_links)
sankey <- sankeyNetwork(
  Links = as.data.frame(links_plot),
  Nodes = nodes,
  Source = "source_id",
  Target = "target_id",
  Value = "value",
  NodeID = "name",
  fontFamily = "Arial",
  fontSize = 12,
  nodeWidth = 24,
  nodePadding = 12,
  sinksRight = FALSE
)
saveWidget(sankey, opt$out_html, selfcontained = TRUE)
message("Wrote ", opt$out_html, " and ", opt$out_links)
