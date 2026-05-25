library(metafor)

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

input_path <- file.path(project_dir, "limma_per_cohort_output", "per_cohort_limma_results.csv")
output_dir <- file.path(project_dir, "exercise_training_meta_output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

per_cohort_all <- read.csv(input_path, check.names = FALSE)

required_cols <- c("gene", "study", "comparison", "logFC", "vi")
missing_cols <- setdiff(required_cols, colnames(per_cohort_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

included_comparisons <- data.frame(
  study = c("GSE116801", "GSE208032", "GSE43471", "GSE58559", "GSE205891"),
  comparison = c(
    "exercise_vs_sedentary",
    "post_training_vs_baseline_paired_rested",
    "exercise_change_vs_control_change",
    "post_intervention_vs_pre_intervention_paired",
    "intensive_lifestyle_change_vs_standard_care_change"
  ),
  inclusion_role = c(
    "exercise_related_cross_sectional",
    "paired_training_response",
    "training_change_vs_control_change",
    "paired_training_plus_diet_response",
    "intensive_lifestyle_change_vs_standard_care_change"
  ),
  note = c(
    "Included as exercise-related, but interpret cautiously because it is not a paired training intervention.",
    "Included as a paired endurance/training before-after response.",
    "Included as training-induced change relative to control change.",
    "Included as a paired post-intervention versus pre-intervention response. GEO/OmicsDI describes exercise training with modest energy deficit; local sample metadata distinguishes diet arms but has no diet-only versus exercise-plus-diet flag, so no GSE58559 samples are excluded.",
    "Included as intensive lifestyle therapy change versus standard care change in subcutaneous adipose tissue."
  ),
  stringsAsFactors = FALSE
)

excluded_comparisons <- data.frame(
  study = character(),
  comparison = character(),
  inclusion_role = character(),
  note = character(),
  stringsAsFactors = FALSE
)

write.csv(
  rbind(included_comparisons, excluded_comparisons),
  file.path(output_dir, "study_inclusion_notes.csv"),
  row.names = FALSE
)

include_key <- paste(included_comparisons$study, included_comparisons$comparison, sep = "__")
row_key <- paste(per_cohort_all$study, per_cohort_all$comparison, sep = "__")

per_cohort <- per_cohort_all[row_key %in% include_key, , drop = FALSE]
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
      all_training_studies = n_studies == n_total_studies,
      method = NA_character_,
      meta_logFC = NA_real_,
      meta_se = NA_real_,
      z = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
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

  data.frame(
    gene = gene_df$gene[[1]],
    n_studies = n_studies,
    all_training_studies = n_studies == n_total_studies,
    method = fit_method,
    meta_logFC = as.numeric(fit$b),
    meta_se = fit$se,
    z = fit$zval,
    p_value = fit$pval,
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
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
gene_results <- gene_results[!vapply(gene_results, is.null, logical(1))]

if (length(gene_results) == 0) {
  stop("No genes were represented in at least two included exercise/training studies.")
}

meta_results <- do.call(rbind, gene_results)
meta_results$adj_p_value <- p.adjust(meta_results$p_value, method = "BH")
meta_results <- meta_results[order(meta_results$p_value), ]

write.csv(
  per_cohort,
  file.path(output_dir, "per_study_differential_expression_exercise_training_only.csv"),
  row.names = FALSE
)

write.csv(
  meta_results,
  file.path(output_dir, "exercise_training_meta_results.csv"),
  row.names = FALSE
)

write.csv(
  meta_results[meta_results$all_training_studies, ],
  file.path(output_dir, "exercise_training_meta_results_all_4_training_related_studies.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "input_rows_after_filter",
    "included_studies",
    "excluded_studies",
    "unique_input_genes",
    "genes_meta_analyzed_at_least_2_training_studies",
    "genes_meta_analyzed_all_5_training_related_studies",
    "fit_errors"
  ),
  value = c(
    nrow(per_cohort),
    n_total_studies,
    length(unique(excluded_comparisons$study)),
    length(unique(per_cohort$gene)),
    nrow(meta_results),
    sum(meta_results$all_training_studies),
    sum(!is.na(meta_results$error))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "exercise_training_meta_summary.csv"),
  row.names = FALSE
)

message(
  "Wrote exercise/training-only meta-analysis results to ",
  file.path(output_dir, "exercise_training_meta_results.csv")
)
print(summary_table)
print(head(meta_results, 20))
