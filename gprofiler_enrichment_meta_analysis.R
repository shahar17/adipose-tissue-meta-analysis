library(gprofiler2)

input_path <- file.path("random_effects_meta_output", "random_effects_meta_results.csv")
output_dir <- "gprofiler_enrichment_output"

adj_p_threshold <- 0.1
i2_threshold <- 50
abs_effect_threshold <- 0.25

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta_results <- read.csv(input_path, check.names = FALSE)

required_cols <- c("gene", "meta_logFC", "z", "p_value", "adj_p_value", "I2")
missing_cols <- setdiff(required_cols, colnames(meta_results))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta_results <- meta_results[
  is.finite(meta_results$meta_logFC) &
    is.finite(meta_results$z) &
    is.finite(meta_results$p_value) &
    is.finite(meta_results$adj_p_value) &
    is.finite(meta_results$I2) &
    nzchar(meta_results$gene),
]

sig <- meta_results[
  meta_results$adj_p_value < adj_p_threshold &
    meta_results$I2 < i2_threshold &
    abs(meta_results$meta_logFC) > abs_effect_threshold,
]
sig_up <- sig[sig$meta_logFC > 0, ]
sig_down <- sig[sig$meta_logFC < 0, ]

gene_sets <- list(
  all_significant = unique(sig$gene),
  upregulated = unique(sig_up$gene),
  downregulated = unique(sig_down$gene)
)

write_gene_list <- function(genes, name) {
  write.csv(
    data.frame(gene = genes, stringsAsFactors = FALSE),
    file.path(output_dir, paste0(name, "_genes.csv")),
    row.names = FALSE
  )
}

for (name in names(gene_sets)) {
  write_gene_list(gene_sets[[name]], name)
}

run_gost <- function(query, name, ordered_query = FALSE) {
  if (length(query) == 0) {
    warning("Skipping ", name, ": no genes.")
    return(data.frame())
  }

  message("Running gProfiler enrichment for ", name, " (", length(query), " genes) ...")

  result <- gost(
    query = query,
    organism = "hsapiens",
    ordered_query = ordered_query,
    significant = TRUE,
    correction_method = "g_SCS",
    sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"),
    evcodes = TRUE
  )

  if (is.null(result) || is.null(result$result) || nrow(result$result) == 0) {
    warning("No significant enrichment returned for ", name)
    return(data.frame())
  }

  out <- result$result
  out$query_name <- name
  out
}

flatten_for_csv <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }

  as.data.frame(lapply(df, function(col) {
    if (is.list(col)) {
      vapply(col, function(x) paste(unlist(x), collapse = ";"), character(1))
    } else {
      col
    }
  }), stringsAsFactors = FALSE, check.names = FALSE)
}

ora_results <- do.call(rbind, lapply(names(gene_sets), function(name) {
  run_gost(gene_sets[[name]], name = name, ordered_query = FALSE)
}))

ora_results <- flatten_for_csv(ora_results)

write.csv(
  ora_results,
  file.path(output_dir, "gprofiler_ora_results.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  analysis = c(
    "adjusted_p_threshold",
    "i2_threshold",
    "absolute_effect_threshold",
    "all_significant_genes",
    "upregulated_significant_genes",
    "downregulated_significant_genes",
    "ora_enriched_terms"
  ),
  value = c(
    adj_p_threshold,
    i2_threshold,
    abs_effect_threshold,
    length(gene_sets$all_significant),
    length(gene_sets$upregulated),
    length(gene_sets$downregulated),
    nrow(ora_results)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "gprofiler_enrichment_summary.csv"),
  row.names = FALSE
)

message("Wrote gProfiler enrichment results to ", output_dir)
print(summary_table)
print(head(ora_results[order(ora_results$p_value), ], 20))
