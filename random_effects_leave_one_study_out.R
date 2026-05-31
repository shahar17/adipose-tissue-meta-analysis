library(metafor)

input_path <- file.path("limma_per_cohort_output", "per_cohort_limma_results.csv")
output_dir <- "leave_one_study_out_meta_output"

adj_p_threshold <- 0.1
i2_threshold <- 50
abs_effect_threshold <- 0.25

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

per_cohort_all <- read.csv(input_path, check.names = FALSE)

required_cols <- c("gene", "study", "comparison", "logFC", "vi")
missing_cols <- setdiff(required_cols, colnames(per_cohort_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

per_cohort_all <- per_cohort_all[
  is.finite(per_cohort_all$logFC) &
    is.finite(per_cohort_all$vi) &
    per_cohort_all$vi > 0,
]

fit_meta_analysis <- function(per_cohort, analysis_name, omitted_study = NA_character_) {
  included_studies <- sort(unique(per_cohort$study))
  n_total_studies <- length(included_studies)

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
        analysis = analysis_name,
        omitted_study = omitted_study,
        gene = gene_df$gene[[1]],
        n_studies = n_studies,
        all_included_studies = n_studies == n_total_studies,
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
      analysis = analysis_name,
      omitted_study = omitted_study,
      gene = gene_df$gene[[1]],
      n_studies = n_studies,
      all_included_studies = n_studies == n_total_studies,
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
    stop("No genes were represented in at least two studies for ", analysis_name)
  }

  meta_results <- do.call(rbind, gene_results)
  meta_results$adj_p_value <- p.adjust(meta_results$p_value, method = "BH")
  meta_results <- meta_results[order(meta_results$p_value), ]

  high_confidence <- meta_results[
    meta_results$adj_p_value < adj_p_threshold &
      meta_results$I2 < i2_threshold &
      abs(meta_results$meta_logFC) > abs_effect_threshold,
  ]

  write.csv(
    meta_results,
    file.path(output_dir, paste0(analysis_name, "_meta_results.csv")),
    row.names = FALSE
  )

  write.csv(
    high_confidence,
    file.path(output_dir, paste0(analysis_name, "_high_confidence_genes.csv")),
    row.names = FALSE
  )

  data.frame(
    analysis = analysis_name,
    omitted_study = omitted_study,
    included_studies = paste(included_studies, collapse = ";"),
    n_included_studies = n_total_studies,
    input_rows = nrow(per_cohort),
    unique_input_genes = length(unique(per_cohort$gene)),
    genes_meta_analyzed_at_least_2_studies = nrow(meta_results),
    genes_meta_analyzed_all_included_studies = sum(meta_results$all_included_studies),
    genes_adj_p_lt_0_1 = sum(meta_results$adj_p_value < adj_p_threshold, na.rm = TRUE),
    high_confidence_genes = nrow(high_confidence),
    high_confidence_up = sum(high_confidence$meta_logFC > 0, na.rm = TRUE),
    high_confidence_down = sum(high_confidence$meta_logFC < 0, na.rm = TRUE),
    fit_errors = sum(!is.na(meta_results$error)),
    stringsAsFactors = FALSE
  )
}

studies <- sort(unique(per_cohort_all$study))

summary_rows <- list(
  fit_meta_analysis(
    per_cohort_all,
    analysis_name = "all_studies",
    omitted_study = NA_character_
  )
)

for (study in studies) {
  per_cohort_subset <- per_cohort_all[per_cohort_all$study != study, , drop = FALSE]
  analysis_name <- paste0("leave_out_", study)

  summary_rows[[analysis_name]] <- fit_meta_analysis(
    per_cohort_subset,
    analysis_name = analysis_name,
    omitted_study = study
  )
}

summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  file.path(output_dir, "leave_one_study_out_summary.csv"),
  row.names = FALSE
)

message("Wrote leave-one-study-out meta-analysis results to ", output_dir)
print(summary_table)
