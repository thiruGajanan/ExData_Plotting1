#' Run ANCOM-BC differential abundance testing
#'
#' @param physeq Phyloseq object.
#' @param formula Model formula.
#' @param p_adj_method Adjustment method.
#' @return Tidy tibble with ANCOM-BC results.
run_ancombc <- function(physeq,
                        formula,
                        p_adj_method = "BH",
                        ...) {
  result <- ANCOMBC::ancombc(
    phyloseq = physeq,
    formula = formula,
    p_adj_method = p_adj_method,
    ...
  )
  primary <- result$res
  beta <- primary$beta |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "beta")
  se <- primary$se |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "se")
  w <- primary$W |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "w")
  p <- primary$p_val |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "p_val")
  q <- primary$q_val |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "q_val")
  diff <- primary$diff_abn |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "coefficient", values_to = "diff_abn")
  beta |>
    dplyr::left_join(se, by = c("feature", "coefficient")) |>
    dplyr::left_join(w, by = c("feature", "coefficient")) |>
    dplyr::left_join(p, by = c("feature", "coefficient")) |>
    dplyr::left_join(q, by = c("feature", "coefficient")) |>
    dplyr::left_join(diff, by = c("feature", "coefficient"))
}

#' Run MaAsLin2 differential abundance modelling
#'
#' @param physeq Phyloseq object.
#' @param metadata Sample metadata tibble.
#' @param fixed_effects Character vector of fixed effects.
#' @param random_effects Optional random effects.
#' @param output_dir Directory for MaAsLin2 outputs.
#' @return Tibble summarising MaAsLin2 results.
run_maaslin2 <- function(physeq,
                         metadata,
                         fixed_effects,
                         random_effects = NULL,
                         output_dir = file.path(tempdir(), "maaslin2"),
                         normalization = "TMM",
                         transform = "LOG") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  otu <- phyloseq::otu_table(physeq)
  if (phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  abund <- as.data.frame(otu)
  metadata_df <- metadata |>
    dplyr::mutate(sample_id = as.character(sample_id)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    tibble::column_to_rownames("sample_id")
  maaslin_fit <- Maaslin2::Maaslin2(
    input_data = abund,
    input_metadata = metadata_df,
    output = output_dir,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    normalization = normalization,
    transform = transform
  )
  readr::read_tsv(file.path(output_dir, "all_results.tsv"), show_col_types = FALSE)
}

#' Summarise the core microbiome
#'
#' @param physeq Phyloseq object.
#' @param detection Minimum relative abundance threshold.
#' @param prevalence Minimum prevalence threshold.
#' @return Tibble listing core taxa and prevalence statistics.
summarize_core_microbiome <- function(physeq,
                                      detection = 0.001,
                                      prevalence = 0.5) {
  core <- microbiome::core(physeq, detection = detection, prevalence = prevalence)
  prev <- microbiome::prevalence(physeq, detection = detection)
  tibble::tibble(
    taxa = names(prev),
    prevalence = as.numeric(prev)
  ) |>
    dplyr::filter(taxa %in% core)
}
