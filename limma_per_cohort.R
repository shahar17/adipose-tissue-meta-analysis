library(limma)

input_dir <- "geo_gene_level_output"
output_dir <- "limma_per_cohort_output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

preprocess_expression <- function(exprs_mat, geo_id) {
  storage.mode(exprs_mat) <- "numeric"

  q99 <- as.numeric(quantile(exprs_mat, 0.99, na.rm = TRUE))

  if (is.finite(q99) && q99 > 50) {
    min_value <- min(exprs_mat, na.rm = TRUE)
    offset <- if (min_value <= 0) abs(min_value) + 1 else 0

    message(
      geo_id,
      ": raw-looking expression scale detected; applying log2(x + ",
      round(offset, 4),
      ") and quantile normalization before limma."
    )

    exprs_mat <- log2(exprs_mat + offset)
    exprs_mat <- normalizeBetweenArrays(exprs_mat, method = "quantile")
  } else {
    message(geo_id, ": expression scale looks log-like; keeping values as downloaded.")
  }

  exprs_mat
}

read_study <- function(geo_id, base_dir = input_dir) {
  study_dir <- file.path(base_dir, geo_id)

  exprs_mat <- as.matrix(
    read.csv(
      file.path(study_dir, "expression_matrix_gene_level.csv"),
      row.names = 1,
      check.names = FALSE
    )
  )
  pheno <- read.csv(
    file.path(study_dir, "sample_metadata.csv"),
    row.names = 1,
    check.names = FALSE
  )

  list(exprs = preprocess_expression(exprs_mat, geo_id), pheno = pheno)
}

extract_participant_id <- function(x) {
  out <- sub("^.*participant_id: *([0-9]+).*$", "\\1", x, ignore.case = TRUE)
  out[!grepl("^[0-9]+$", out)] <- NA_character_
  out
}

extract_regex_group <- function(x, pattern) {
  matches <- regexec(pattern, x, ignore.case = TRUE)
  out <- regmatches(x, matches)

  vapply(out, function(hit) {
    if (length(hit) >= 2) {
      hit[[2]]
    } else {
      NA_character_
    }
  }, character(1))
}

clean_limma_table <- function(tt, geo_id, comparison, n_case, n_control) {
  se <- abs(tt$logFC / tt$t)

  data.frame(
    gene = rownames(tt),
    study = geo_id,
    comparison = comparison,
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    p_value = tt$P.Value,
    adj_p_value = tt$adj.P.Val,
    B = tt$B,
    se = se,
    vi = se^2,
    n_case = n_case,
    n_control = n_control,
    stringsAsFactors = FALSE
  )
}

run_paired_change_limma <- function(exprs_mat, pheno, geo_id, subject_col, timepoint_col,
                                    baseline_label, post_label, comparison) {
  keep <- pheno[[timepoint_col]] %in% c(baseline_label, post_label) &
    !is.na(pheno[[subject_col]])

  exprs_mat <- exprs_mat[, keep, drop = FALSE]
  pheno <- pheno[keep, , drop = FALSE]

  sample_key <- paste(pheno[[subject_col]], pheno[[timepoint_col]], sep = "__")
  averaged <- t(rowsum(t(exprs_mat), group = sample_key, reorder = FALSE))
  key_parts <- strsplit(colnames(averaged), "__", fixed = TRUE)

  averaged_pheno <- data.frame(
    sample_key = colnames(averaged),
    subject_id = vapply(key_parts, `[`, character(1), 1),
    timepoint = vapply(key_parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  by_subject <- split(averaged_pheno, averaged_pheno$subject_id)
  change_list <- list()

  for (subject_id in names(by_subject)) {
    rows <- by_subject[[subject_id]]
    baseline_key <- rows$sample_key[rows$timepoint == baseline_label]
    post_key <- rows$sample_key[rows$timepoint == post_label]

    if (length(baseline_key) == 1 && length(post_key) == 1) {
      change_list[[subject_id]] <- averaged[, post_key] - averaged[, baseline_key]
    }
  }

  change_mat <- do.call(cbind, change_list)
  colnames(change_mat) <- names(change_list)

  design <- matrix(1, nrow = ncol(change_mat), ncol = 1)
  colnames(design) <- "post_minus_baseline"

  fit <- lmFit(change_mat, design)
  fit <- eBayes(fit)

  tt <- topTable(fit, coef = "post_minus_baseline", number = Inf, sort.by = "none")
  clean_limma_table(
    tt,
    geo_id = geo_id,
    comparison = comparison,
    n_case = ncol(change_mat),
    n_control = NA_integer_
  )
}

run_gse43471 <- function() {
  geo_id <- "GSE43471"
  dat <- read_study(geo_id)

  pheno <- dat$pheno[colnames(dat$exprs), , drop = FALSE]
  pheno$group <- pheno[["sample group:ch1"]]
  pheno$participant_id <- extract_participant_id(pheno$characteristics_ch1)
  pheno$timepoint <- ifelse(
    grepl("baseline", pheno$source_name_ch1, ignore.case = TRUE),
    "baseline",
    ifelse(grepl("6 month", pheno$source_name_ch1, ignore.case = TRUE), "post_training", NA)
  )

  keep <- pheno$group %in% c("control", "exercise") &
    pheno$timepoint %in% c("baseline", "post_training") &
    !is.na(pheno$participant_id)

  exprs_mat <- dat$exprs[, keep, drop = FALSE]
  pheno <- pheno[keep, , drop = FALSE]

  sample_key <- paste(pheno$group, pheno$participant_id, pheno$timepoint, sep = "__")
  averaged <- t(rowsum(t(exprs_mat), group = sample_key, reorder = FALSE))
  key_parts <- strsplit(colnames(averaged), "__", fixed = TRUE)

  averaged_pheno <- data.frame(
    sample_key = colnames(averaged),
    group = vapply(key_parts, `[`, character(1), 1),
    participant_id = vapply(key_parts, `[`, character(1), 2),
    timepoint = vapply(key_parts, `[`, character(1), 3),
    stringsAsFactors = FALSE
  )

  by_subject <- split(averaged_pheno, averaged_pheno$participant_id)
  change_list <- list()
  change_group <- character()

  for (subject_id in names(by_subject)) {
    rows <- by_subject[[subject_id]]
    baseline_key <- rows$sample_key[rows$timepoint == "baseline"]
    post_key <- rows$sample_key[rows$timepoint == "post_training"]

    if (length(baseline_key) == 1 && length(post_key) == 1) {
      change_list[[subject_id]] <- averaged[, post_key] - averaged[, baseline_key]
      change_group <- c(change_group, rows$group[1])
    }
  }

  change_mat <- do.call(cbind, change_list)
  colnames(change_mat) <- names(change_list)
  group <- factor(change_group, levels = c("control", "exercise"))

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrast <- makeContrasts(exercise - control, levels = design)

  fit <- lmFit(change_mat, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)

  tt <- topTable(fit, number = Inf, sort.by = "none")
  clean_limma_table(
    tt,
    geo_id = geo_id,
    comparison = "exercise_change_vs_control_change",
    n_case = sum(group == "exercise"),
    n_control = sum(group == "control")
  )
}

run_gse208032 <- function() {
  geo_id <- "GSE208032"
  dat <- read_study(geo_id)

  pheno <- dat$pheno[colnames(dat$exprs), , drop = FALSE]
  pheno$subject_id <- extract_regex_group(pheno$source_name_ch1, "subject ([0-9]+)")
  pheno$timepoint <- ifelse(
    grepl("untrained before exercise", pheno$source_name_ch1, ignore.case = TRUE),
    "baseline",
    ifelse(
      grepl("post training before exercise", pheno$source_name_ch1, ignore.case = TRUE),
      "post_training",
      NA_character_
    )
  )

  run_paired_change_limma(
    exprs_mat = dat$exprs,
    pheno = pheno,
    geo_id = geo_id,
    subject_col = "subject_id",
    timepoint_col = "timepoint",
    baseline_label = "baseline",
    post_label = "post_training",
    comparison = "post_training_vs_baseline_paired_rested"
  )
}

run_gse58559 <- function() {
  geo_id <- "GSE58559"
  dat <- read_study(geo_id)

  pheno <- dat$pheno[colnames(dat$exprs), , drop = FALSE]
  pheno$subject_id <- pheno[["individual:ch1"]]
  pheno$timepoint <- pheno[["time point:ch1"]]

  run_paired_change_limma(
    exprs_mat = dat$exprs,
    pheno = pheno,
    geo_id = geo_id,
    subject_col = "subject_id",
    timepoint_col = "timepoint",
    baseline_label = "pre-diet",
    post_label = "post-diet",
    comparison = "post_intervention_vs_pre_intervention_paired"
  )
}

run_gse116801 <- function() {
  geo_id <- "GSE116801"
  dat <- read_study(geo_id)

  pheno <- dat$pheno[colnames(dat$exprs), , drop = FALSE]
  group <- factor(tolower(pheno$source_name_ch1), levels = c("sedentary", "exercise"))

  keep <- !is.na(group)
  exprs_mat <- dat$exprs[, keep, drop = FALSE]
  group <- droplevels(group[keep])

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrast <- makeContrasts(exercise - sedentary, levels = design)

  fit <- lmFit(exprs_mat, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)

  tt <- topTable(fit, number = Inf, sort.by = "none")
  clean_limma_table(
    tt,
    geo_id = geo_id,
    comparison = "exercise_vs_sedentary",
    n_case = sum(group == "exercise"),
    n_control = sum(group == "sedentary")
  )
}

study_notes <- data.frame(
  geo_id = c("GSE43471", "GSE116801", "GSE58559", "GSE208032", "GSE205891"),
  status = c("included", "included", "included", "included", "excluded"),
  reason = c(
    "Mapped to gene level and has control plus exercise samples over time.",
    "Mapped to gene level and has exercise versus sedentary groups.",
    "Mapped to gene level and included as paired post-intervention versus pre-intervention comparison. GEO/OmicsDI describes exercise training with modest energy deficit; local sample metadata distinguishes diet arms but has no diet-only versus exercise-plus-diet flag.",
    "Mapped through GPL23159 SPOT_ID.1 annotations and included as paired rested post-training versus baseline comparison.",
    "Series matrix contained zero expression-feature rows."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  study_notes,
  file.path(output_dir, "study_analysis_notes.csv"),
  row.names = FALSE
)

results <- rbind(
  run_gse43471(),
  run_gse116801(),
  run_gse58559(),
  run_gse208032()
)

write.csv(
  results,
  file.path(output_dir, "per_cohort_limma_results.csv"),
  row.names = FALSE
)

split_results <- split(results, results$study)
for (study_id in names(split_results)) {
  write.csv(
    split_results[[study_id]],
    file.path(output_dir, paste0(study_id, "_limma_results.csv")),
    row.names = FALSE
  )
}

message("Wrote limma results to ", file.path(output_dir, "per_cohort_limma_results.csv"))
print(head(results, 20))
