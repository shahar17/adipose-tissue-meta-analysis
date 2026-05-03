script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

input_dir <- file.path(project_dir, "geo_extract_output")
output_dir <- file.path(project_dir, "geo_gene_level_output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pick_mapping_column <- function(feature_df) {
  candidates <- c(
    "Gene Symbol",
    "Gene symbol",
    "Symbol",
    "Gene.Symbol",
    "GENE_SYMBOL",
    "ENTREZ_GENE_ID",
    "Entrez_Gene_ID",
    "RefSeq Transcript ID",
    "RefSeq_ID",
    "SPOT_ID.1",
    "SPOT_ID"
  )

  hits <- intersect(candidates, colnames(feature_df))
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- sub(" ///.*$", "", x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x[x == "---"] <- NA_character_
  x
}

extract_parenthetical_symbol <- function(x) {
  x <- as.character(x)
  matches <- gregexpr("\\(([A-Z][A-Z0-9.-]{1,})\\)", x, perl = TRUE)
  out <- regmatches(x, matches)

  vapply(out, function(hit) {
    if (length(hit) == 0 || identical(hit, character(0))) {
      return(NA_character_)
    }

    symbols <- gsub("^\\(|\\)$", "", hit)
    symbols <- symbols[!symbols %in% c("RNA", "DNA", "HGNC", "SOURCE")]
    if (length(symbols) == 0) {
      NA_character_
    } else {
      symbols[[1]]
    }
  }, character(1))
}

derive_gene_ids <- function(feature_df, mapping_col) {
  if (mapping_col %in% c("SPOT_ID.1", "SPOT_ID")) {
    return(clean_gene_id(extract_parenthetical_symbol(feature_df[[mapping_col]])))
  }

  clean_gene_id(feature_df[[mapping_col]])
}

read_feature_table <- function(path) {
  out <- tryCatch(
    read.csv(path, row.names = 1, check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(data.frame())
  }

  out
}

collapse_to_gene_level <- function(exprs_mat, gene_ids) {
  keep <- !is.na(gene_ids)
  exprs_mat <- exprs_mat[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]

  collapsed_sum <- rowsum(exprs_mat, group = gene_ids, reorder = FALSE)
  counts <- as.vector(table(factor(gene_ids, levels = rownames(collapsed_sum))))
  collapsed_sum / counts
}

map_one_study <- function(geo_id, base_input = input_dir, base_output = output_dir) {
  study_dir <- file.path(base_input, geo_id)
  exprs_path <- file.path(study_dir, "expression_matrix.csv")
  feature_path <- file.path(study_dir, "feature_metadata.csv")
  sample_path <- file.path(study_dir, "sample_metadata.csv")

  exprs_mat <- read.csv(exprs_path, row.names = 1, check.names = FALSE)
  feature_df <- read_feature_table(feature_path)
  sample_df <- read.csv(sample_path, row.names = 1, check.names = FALSE)

  if (nrow(feature_df) == 0) {
    return(data.frame(
      geo_id = geo_id,
      mapping_status = "empty_feature_metadata",
      chosen_mapping_column = NA_character_,
      n_input_features = nrow(exprs_mat),
      n_mapped_features = NA_integer_,
      n_gene_rows = NA_integer_,
      gene_level_expression_path = NA_character_,
      sample_metadata_path = sample_path,
      error = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  mapping_col <- pick_mapping_column(feature_df)

  if (is.na(mapping_col)) {
    return(data.frame(
      geo_id = geo_id,
      mapping_status = "no_mapping_column_found",
      chosen_mapping_column = NA_character_,
      n_input_features = nrow(exprs_mat),
      n_mapped_features = NA_integer_,
      n_gene_rows = NA_integer_,
      gene_level_expression_path = NA_character_,
      sample_metadata_path = sample_path,
      error = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  common_ids <- intersect(rownames(exprs_mat), rownames(feature_df))
  exprs_mat <- as.matrix(exprs_mat[common_ids, , drop = FALSE])
  feature_df <- feature_df[common_ids, , drop = FALSE]

  gene_ids <- derive_gene_ids(feature_df, mapping_col)
  gene_level <- collapse_to_gene_level(exprs_mat, gene_ids)

  out_dir <- file.path(base_output, geo_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    gene_level,
    file.path(out_dir, "expression_matrix_gene_level.csv"),
    row.names = TRUE
  )
  write.csv(
    sample_df,
    file.path(out_dir, "sample_metadata.csv"),
    row.names = TRUE
  )

  data.frame(
    geo_id = geo_id,
    mapping_status = "mapped",
    chosen_mapping_column = mapping_col,
    n_input_features = length(common_ids),
    n_mapped_features = sum(!is.na(gene_ids)),
    n_gene_rows = nrow(gene_level),
    gene_level_expression_path = file.path(out_dir, "expression_matrix_gene_level.csv"),
    sample_metadata_path = file.path(out_dir, "sample_metadata.csv"),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_mapping <- function(base_input = input_dir, base_output = output_dir) {
  geo_ids <- list.dirs(base_input, recursive = FALSE, full.names = FALSE)
  geo_ids <- geo_ids[nzchar(geo_ids)]

  results <- lapply(geo_ids, function(geo_id) {
    tryCatch(
      map_one_study(geo_id, base_input = base_input, base_output = base_output),
      error = function(e) {
        data.frame(
          geo_id = geo_id,
          mapping_status = "error",
          chosen_mapping_column = NA_character_,
          n_input_features = NA_integer_,
          n_mapped_features = NA_integer_,
          n_gene_rows = NA_integer_,
          gene_level_expression_path = NA_character_,
          sample_metadata_path = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  })

  manifest <- do.call(rbind, results)

  write.csv(
    manifest,
    file.path(base_output, "gene_mapping_manifest.csv"),
    row.names = FALSE
  )

  message("Finished. Manifest written to ", file.path(base_output, "gene_mapping_manifest.csv"))
  manifest
}

manifest <- run_mapping()
print(manifest)
