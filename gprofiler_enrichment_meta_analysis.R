library(gprofiler2)

input_path <- file.path("random_effects_meta_output", "random_effects_meta_results.csv")
output_dir <- "gprofiler_enrichment_output"

strict_adj_p_threshold <- 0.1
strict_i2_threshold <- 50
strict_abs_effect_threshold <- 0.25

relaxed_adj_p_threshold <- 0.1
relaxed_abs_effect_threshold <- 0.1

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

make_gene_sets <- function(sig) {
  sig_up <- sig[sig$meta_logFC > 0, ]
  sig_down <- sig[sig$meta_logFC < 0, ]

  list(
    all_significant = unique(sig$gene),
    upregulated = unique(sig_up$gene),
    downregulated = unique(sig_down$gene)
  )
}

strict_sig <- meta_results[
  meta_results$adj_p_value < strict_adj_p_threshold &
    meta_results$I2 < strict_i2_threshold &
    abs(meta_results$meta_logFC) > strict_abs_effect_threshold,
]

relaxed_sig <- meta_results[
  meta_results$adj_p_value < relaxed_adj_p_threshold &
    abs(meta_results$meta_logFC) > relaxed_abs_effect_threshold,
]

strict_gene_sets <- make_gene_sets(strict_sig)
relaxed_gene_sets <- make_gene_sets(relaxed_sig)

write_gene_list <- function(genes, name) {
  write.csv(
    data.frame(gene = genes, stringsAsFactors = FALSE),
    file.path(output_dir, paste0(name, "_genes.csv")),
    row.names = FALSE
  )
}

write_gene_set_group <- function(gene_sets, prefix) {
  for (name in names(gene_sets)) {
    write_gene_list(gene_sets[[name]], paste0(prefix, "_", name))
  }
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

run_gene_set_group <- function(gene_sets, prefix) {
  results <- do.call(rbind, lapply(names(gene_sets), function(name) {
    run_gost(gene_sets[[name]], name = paste0(prefix, "_", name), ordered_query = FALSE)
  }))

  flatten_for_csv(results)
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

write_gene_set_group(strict_gene_sets, "strict")
write_gene_set_group(relaxed_gene_sets, "relaxed_no_i2_abs0.1")

for (name in names(strict_gene_sets)) {
  write_gene_list(strict_gene_sets[[name]], name)
}

strict_ora_results <- run_gene_set_group(strict_gene_sets, "strict")
relaxed_ora_results <- run_gene_set_group(relaxed_gene_sets, "relaxed_no_i2_abs0.1")

ora_results <- rbind(strict_ora_results, relaxed_ora_results)

write.csv(
  strict_ora_results,
  file.path(output_dir, "gprofiler_ora_results_strict.csv"),
  row.names = FALSE
)

write.csv(
  relaxed_ora_results,
  file.path(output_dir, "gprofiler_ora_results_relaxed_no_i2_abs0.1.csv"),
  row.names = FALSE
)

write.csv(
  strict_ora_results,
  file.path(output_dir, "gprofiler_ora_results.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  analysis = c(
    "strict_adjusted_p_threshold",
    "strict_i2_threshold",
    "strict_absolute_effect_threshold",
    "strict_all_significant_genes",
    "strict_upregulated_significant_genes",
    "strict_downregulated_significant_genes",
    "strict_ora_enriched_terms",
    "relaxed_adjusted_p_threshold",
    "relaxed_i2_threshold",
    "relaxed_absolute_effect_threshold",
    "relaxed_all_significant_genes",
    "relaxed_upregulated_significant_genes",
    "relaxed_downregulated_significant_genes",
    "relaxed_ora_enriched_terms"
  ),
  value = c(
    strict_adj_p_threshold,
    strict_i2_threshold,
    strict_abs_effect_threshold,
    length(strict_gene_sets$all_significant),
    length(strict_gene_sets$upregulated),
    length(strict_gene_sets$downregulated),
    nrow(strict_ora_results),
    relaxed_adj_p_threshold,
    "not_applied",
    relaxed_abs_effect_threshold,
    length(relaxed_gene_sets$all_significant),
    length(relaxed_gene_sets$upregulated),
    length(relaxed_gene_sets$downregulated),
    nrow(relaxed_ora_results)
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
