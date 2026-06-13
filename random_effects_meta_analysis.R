library(metafor)

input_path <- file.path("limma_per_cohort_output", "per_cohort_limma_results.csv")
output_dir <- "random_effects_meta_output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

per_cohort <- read.csv(input_path, check.names = FALSE)

required_cols <- c("gene", "study", "comparison", "logFC", "vi")
missing_cols <- setdiff(required_cols, colnames(per_cohort))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

per_cohort <- per_cohort[
  is.finite(per_cohort$logFC) &
    is.finite(per_cohort$vi) &
    per_cohort$vi > 0,
]

all_studies <- sort(unique(per_cohort$study))
n_total_studies <- length(all_studies)

fit_gene <- function(gene_df) {
  gene_df <- gene_df[order(gene_df$study), , drop = FALSE]
  n_studies <- length(unique(gene_df$study))

  if (n_studies < 2) {
    return(NULL)
  }

  fit_method <- "REML"
  fit <- tryCatch(
    rma.uni(
      yi = gene_df$logFC,
      vi = gene_df$vi,
      method = fit_method,
      test = "z"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    fit_method <- "DL"
    fit <- tryCatch(
      rma.uni(
        yi = gene_df$logFC,
        vi = gene_df$vi,
        method = fit_method,
        test = "z"
      ),
      error = function(e) e
    )
  }

  if (inherits(fit, "error")) {
    return(data.frame(
      gene = gene_df$gene[[1]],
      n_studies = n_studies,
      all_studies = n_studies == n_total_studies,
      method = NA_character_,
      meta_logFC = NA_real_,
      meta_se = NA_real_,
      z = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      prediction_lower = NA_real_,
      prediction_upper = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      H2 = NA_real_,
      Q = NA_real_,
      Q_p_value = NA_real_,
      studies = paste(unique(gene_df$study), collapse = ";"),
      comparisons = paste(unique(gene_df$comparison), collapse = ";"),
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  predicted <- predict(fit)

  data.frame(
    gene = gene_df$gene[[1]],
    n_studies = n_studies,
    all_studies = n_studies == n_total_studies,
    method = fit_method,
    meta_logFC = as.numeric(fit$b),
    meta_se = fit$se,
    z = fit$zval,
    p_value = fit$pval,
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
    prediction_lower = predicted$pi.lb,
    prediction_upper = predicted$pi.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    H2 = fit$H2,
    Q = fit$QE,
    Q_p_value = fit$QEp,
    studies = paste(unique(gene_df$study), collapse = ";"),
    comparisons = paste(unique(gene_df$comparison), collapse = ";"),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

gene_results <- lapply(split(per_cohort, per_cohort$gene), fit_gene)
meta_results <- do.call(rbind, gene_results[!vapply(gene_results, is.null, logical(1))])

meta_results$adj_p_value <- p.adjust(meta_results$p_value, method = "BH")
meta_results <- meta_results[order(meta_results$p_value), ]

write.csv(
  meta_results,
  file.path(output_dir, "random_effects_meta_results.csv"),
  row.names = FALSE
)

write.csv(
  meta_results[meta_results$all_studies, ],
  file.path(output_dir, paste0("random_effects_meta_results_all_", n_total_studies, "_studies.csv")),
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "input_rows",
    "input_studies",
    "unique_input_genes",
    "genes_meta_analyzed_at_least_2_studies",
    paste0("genes_meta_analyzed_all_", n_total_studies, "_studies"),
    "fit_errors"
  ),
  value = c(
    nrow(per_cohort),
    n_total_studies,
    length(unique(per_cohort$gene)),
    nrow(meta_results),
    sum(meta_results$all_studies),
    sum(!is.na(meta_results$error))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "random_effects_meta_summary.csv"),
  row.names = FALSE
)

message(
  "Wrote random-effects meta-analysis results to ",
  file.path(output_dir, "random_effects_meta_results.csv")
)
print(summary_table)
print(head(meta_results, 20))
