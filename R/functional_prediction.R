#' Predict functional profiles from 16S ASVs
#'
#' @param physeq Phyloseq object with ASV counts and taxonomy.
#' @param method Functional prediction backend (`picrust2` or `tax4fun2`).
#' @param output_dir Directory for intermediate and result files.
#' @param picrust2_pipeline Path or command name for `picrust2_pipeline.py`.
#' @param threads Number of threads for PICRUSt2.
#' @param reference_data Optional reference dataset required for Tax4Fun2.
#' @return List containing pathway and KO abundance tables plus metadata.
predict_functional_profiles <- function(physeq,
                                        method = c("picrust2", "tax4fun2"),
                                        output_dir = file.path("outputs", "functional"),
                                        picrust2_pipeline = "picrust2_pipeline.py",
                                        threads = 1,
                                        reference_data = NULL) {
  method <- match.arg(method)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (method == "picrust2") {
    picrust2_exec <- Sys.which(picrust2_pipeline)
    if (picrust2_exec == "") {
      warning("PICRUSt2 pipeline not found. Skipping functional prediction.")
      return(list(
        method = method,
        pathway_abundance = tibble::tibble(),
        ko_abundance = tibble::tibble(),
        output_dir = output_dir
      ))
    }
    otu <- phyloseq::otu_table(physeq)
    if (!phyloseq::taxa_are_rows(physeq)) {
      otu <- t(otu)
    }
    seqs <- Biostrings::DNAStringSet(rownames(otu))
    fasta_path <- file.path(output_dir, "asv_sequences.fasta")
    Biostrings::writeXStringSet(seqs, fasta_path)
    biom_path <- file.path(output_dir, "asv_table.biom")
    biom <- biomformat::make_biom(data = otu)
    biomformat::write_biom(biom, biom_path)
    args <- c(
      "-s", fasta_path,
      "-i", biom_path,
      "-o", output_dir,
      "-p", as.character(threads),
      "--stratified"
    )
    message("Running PICRUSt2 pipeline...")
    status <- system2(picrust2_exec, args = args)
    if (!identical(status, 0L)) {
      stop("PICRUSt2 pipeline failed with status ", status)
    }
    pathway_file <- file.path(output_dir, "pathways_out", "path_abun_unstrat.tsv")
    ko_file <- file.path(output_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv")
    pathway_tbl <- if (file.exists(pathway_file)) {
      readr::read_tsv(pathway_file, comment = "#", show_col_types = FALSE)
    } else {
      tibble::tibble()
    }
    ko_tbl <- if (file.exists(ko_file)) {
      readr::read_tsv(ko_file, comment = "#", show_col_types = FALSE)
    } else {
      tibble::tibble()
    }
    list(
      method = method,
      pathway_abundance = pathway_tbl,
      ko_abundance = ko_tbl,
      output_dir = output_dir
    )
  } else {
    if (!requireNamespace("Tax4Fun2", quietly = TRUE)) {
      stop("Tax4Fun2 package is not installed. Install it or use method = 'picrust2'.")
    }
    if (is.null(reference_data)) {
      stop("`reference_data` must be provided for Tax4Fun2 predictions.")
    }
    otu <- phyloseq::otu_table(physeq)
    if (!phyloseq::taxa_are_rows(physeq)) {
      otu <- t(otu)
    }
    tax <- as.data.frame(phyloseq::tax_table(physeq))
    tax$Sequence <- rownames(tax)
    tax_input <- list(
      otu_table = as.matrix(otu),
      tax_table = tax
    )
    prediction <- Tax4Fun2::runFun(
      tax_input,
      referenceData = reference_data,
      shortReadMode = TRUE,
      normCopyNo = TRUE
    )
    list(
      method = method,
      pathway_abundance = tibble::as_tibble(prediction$path_abundance, rownames = "pathway"),
      ko_abundance = tibble::as_tibble(prediction$functional_profiles, rownames = "ko"),
      output_dir = output_dir
    )
  }
}

#' Analyse pathway differences between groups
#'
#' @param pathway_abundance Tibble of predicted pathway abundances.
#' @param metadata Sample metadata tibble.
#' @param group Column in metadata representing group labels.
#' @param top_n Number of top pathways to plot.
#' @return A list containing summary statistics and a ggplot object.
analyze_pathway_differences <- function(pathway_abundance,
                                        metadata,
                                        group,
                                        top_n = 20) {
  if (nrow(pathway_abundance) == 0) {
    warning("Pathway abundance table is empty; skipping pathway analysis.")
    return(list(summary = tibble::tibble(), plot = ggplot2::ggplot()))
  }
  group <- rlang::ensym(group)
  pathway_tbl <- pathway_abundance
  names(pathway_tbl)[1] <- "pathway"
  tidy <- pathway_tbl |>
    tidyr::pivot_longer(-pathway, names_to = "sample_id", values_to = "abundance") |>
    dplyr::left_join(metadata, by = "sample_id") |>
    dplyr::filter(!is.na(!!group))
  summary <- tidy |>
    dplyr::group_by(pathway, !!group) |>
    dplyr::summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(pathway) |>
    dplyr::mutate(range = max(mean_abundance) - min(mean_abundance)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(range)) |>
    dplyr::slice_head(n = top_n * dplyr::n_distinct(!!group))
  plot <- summary |>
    ggplot2::ggplot(ggplot2::aes(x = reorder(pathway, range), y = mean_abundance, fill = !!group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Pathway", y = "Mean predicted abundance", fill = rlang::as_name(group)) +
    ggplot2::theme_minimal()
  list(summary = summary, plot = plot)
}
