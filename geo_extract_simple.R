library(GEOquery)
library(Biobase)

# Edit this list with the GEO series you want to extract.
geo_ids <- c(
  "GSE43471",
  "GSE208032",
  "GSE116801",
  "GSE58559",
  "GSE205891"
)
geo_ids_env <- Sys.getenv("GEO_IDS")
if (nzchar(geo_ids_env)) {
  geo_ids <- trimws(strsplit(geo_ids_env, ",", fixed = TRUE)[[1]])
}

output_dir <- "geo_extract_output"

pick_eset <- function(gset) {
  if (inherits(gset, "ExpressionSet")) {
    return(gset)
  }

  if (length(gset) == 0) {
    stop("No ExpressionSet objects returned.")
  }

  if (length(gset) == 1) {
    return(gset[[1]])
  }

  sizes <- vapply(gset, ncol, integer(1))
  gset[[which.max(sizes)]]
}

safe_path <- function(...) {
  path <- file.path(...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

extract_gse205891_supplement <- function(base_dir = output_dir) {
  geo_id <- "GSE205891"
  study_dir <- file.path(base_dir, geo_id)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

  xlsx_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205891/suppl/GSE205891_LID_all.gene_moderated_log2cpm.xlsx"
  xlsx_path <- file.path(study_dir, "GSE205891_LID_all.gene_moderated_log2cpm.xlsx")
  if (!file.exists(xlsx_path)) {
    message("Downloading ", geo_id, " supplementary XLSX ...")
    utils::download.file(xlsx_url, xlsx_path, mode = "wb", quiet = FALSE)
  }

  message("Reading ", geo_id, " supplementary XLSX ...")
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read ", basename(xlsx_path), ". Install it with install.packages('readxl').")
  }
  workbook <- as.data.frame(readxl::read_excel(xlsx_path), check.names = FALSE)
  gene_symbols <- workbook[["Gene symbol"]]
  fat_cols <- grep("^Fat_", colnames(workbook), value = TRUE)
  if (length(fat_cols) == 0) {
    stop("No Fat_* columns found in ", basename(xlsx_path))
  }

  exprs_mat <- as.matrix(workbook[, fat_cols, drop = FALSE])
  storage.mode(exprs_mat) <- "numeric"
  rownames(exprs_mat) <- make.unique(gene_symbols)

  pheno_path <- file.path(study_dir, "sample_metadata.csv")
  if (file.exists(pheno_path) && file.info(pheno_path)$size > 3) {
    pheno <- read.csv(pheno_path, row.names = 1, check.names = FALSE)
  } else {
    gset_raw <- getGEO(geo_id, GSEMatrix = TRUE)
    eset <- pick_eset(gset_raw)
    pheno <- pData(eset)
  }
  pheno <- pheno[pheno$title %in% fat_cols, , drop = FALSE]
  pheno <- pheno[match(fat_cols, pheno$title), , drop = FALSE]

  colnames(exprs_mat) <- rownames(pheno)

  feature <- data.frame(
    `Gene Symbol` = gene_symbols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(feature) <- rownames(exprs_mat)

  write.csv(
    exprs_mat,
    safe_path(study_dir, "expression_matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    pheno,
    safe_path(study_dir, "sample_metadata.csv"),
    row.names = TRUE
  )
  write.csv(
    feature,
    safe_path(study_dir, "feature_metadata.csv"),
    row.names = TRUE
  )

  data.frame(
    geo_id = geo_id,
    platform_id = "supplementary_gene_moderated_log2cpm_xlsx",
    n_features = nrow(exprs_mat),
    n_samples = ncol(exprs_mat),
    expression_matrix_path = file.path(study_dir, "expression_matrix.csv"),
    sample_metadata_path = file.path(study_dir, "sample_metadata.csv"),
    feature_metadata_path = file.path(study_dir, "feature_metadata.csv"),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

extract_one_study <- function(geo_id, base_dir = output_dir) {
  if (identical(geo_id, "GSE205891")) {
    return(extract_gse205891_supplement(base_dir = base_dir))
  }

  message("Downloading ", geo_id, " ...")

  gset_raw <- getGEO(geo_id, GSEMatrix = TRUE)
  eset <- pick_eset(gset_raw)

  exprs_mat <- exprs(eset)
  pheno <- pData(eset)
  feature <- fData(eset)
  platform_id <- if (length(annotation(eset)) > 0) as.character(annotation(eset)[[1]]) else NA_character_

  study_dir <- file.path(base_dir, geo_id)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    exprs_mat,
    safe_path(study_dir, "expression_matrix.csv"),
    row.names = TRUE
  )
  write.csv(
    pheno,
    safe_path(study_dir, "sample_metadata.csv"),
    row.names = TRUE
  )
  write.csv(
    feature,
    safe_path(study_dir, "feature_metadata.csv"),
    row.names = TRUE
  )

  data.frame(
    geo_id = geo_id,
    platform_id = platform_id,
    n_features = nrow(exprs_mat),
    n_samples = ncol(exprs_mat),
    expression_matrix_path = file.path(study_dir, "expression_matrix.csv"),
    sample_metadata_path = file.path(study_dir, "sample_metadata.csv"),
    feature_metadata_path = file.path(study_dir, "feature_metadata.csv"),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_extract <- function(ids = geo_ids, base_dir = output_dir) {
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(ids, function(geo_id) {
    tryCatch(
      extract_one_study(geo_id, base_dir = base_dir),
      error = function(e) {
        data.frame(
          geo_id = geo_id,
          platform_id = NA_character_,
          n_features = NA_integer_,
          n_samples = NA_integer_,
          expression_matrix_path = NA_character_,
          sample_metadata_path = NA_character_,
          feature_metadata_path = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  })

  manifest <- do.call(rbind, results)

  write.csv(
    manifest,
    file.path(base_dir, "study_manifest.csv"),
    row.names = FALSE
  )

  message("Finished. Manifest written to ", file.path(base_dir, "study_manifest.csv"))
  manifest
}

manifest <- run_extract()
print(manifest)
