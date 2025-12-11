#' Load sample metadata
#'
#' @param metadata_path Path to a tab-delimited metadata file.
#' @return A tibble containing sample metadata with `sample_id` as character.
load_metadata <- function(metadata_path) {
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.x), NA, .x))) |>
    dplyr::mutate(sample_id = as.character(sample_id))
  if (anyDuplicated(metadata$sample_id) > 0) {
    stop("Sample identifiers must be unique in the metadata table.")
  }
  metadata
}

#' Identify paired FASTQ files in a directory
#'
#' @param fastq_dir Directory containing raw FASTQ files.
#' @param forward_suffix Pattern identifying forward reads.
#' @param reverse_suffix Pattern identifying reverse reads.
#' @return A tibble linking sample IDs to forward and reverse FASTQ files.
find_fastqs <- function(fastq_dir,
                        forward_suffix = "_R1_001.fastq.gz",
                        reverse_suffix = "_R2_001.fastq.gz") {
  if (!dir.exists(fastq_dir)) {
    stop("FASTQ directory does not exist: ", fastq_dir)
  }
  forward <- sort(Sys.glob(file.path(fastq_dir, paste0("*", forward_suffix))))
  reverse <- sort(Sys.glob(file.path(fastq_dir, paste0("*", reverse_suffix))))
  if (length(forward) == 0) {
    warning("No forward FASTQ files detected in ", fastq_dir)
  }
  if (length(forward) != length(reverse)) {
    stop("Mismatch in the number of forward and reverse FASTQ files.")
  }
  sample_ids <- basename(forward) |>
    stringr::str_remove(forward_suffix)
  tibble::tibble(
    sample_id = sample_ids,
    forward = forward,
    reverse = reverse
  )
}

#' Generate read quality profiles for a subset of samples
#'
#' @param fastq_table Tibble returned by `find_fastqs`.
#' @param output_dir Directory where quality profile plots will be saved.
#' @param max_samples Maximum number of samples to visualise.
#' @return Tibble with file paths to quality profile images for each sample.
#' @details Plots are generated using `dada2::plotQualityProfile` and saved as PNGs.
generate_quality_profiles <- function(fastq_table,
                                       output_dir = file.path("outputs", "quality_profiles"),
                                       max_samples = 6) {
  if (nrow(fastq_table) == 0) {
    return(tibble::tibble())
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  subset_tbl <- fastq_table |>
    dplyr::slice_head(n = max_samples)
  purrr::pmap_dfr(
    subset_tbl,
    function(sample_id, forward, reverse) {
      forward_plot <- dada2::plotQualityProfile(forward) +
        ggplot2::ggtitle(paste(sample_id, "Forward"))
      reverse_plot <- dada2::plotQualityProfile(reverse) +
        ggplot2::ggtitle(paste(sample_id, "Reverse"))
      forward_path <- file.path(output_dir, paste0(sample_id, "_forward.png"))
      reverse_path <- file.path(output_dir, paste0(sample_id, "_reverse.png"))
      ggplot2::ggsave(forward_path, forward_plot, width = 8, height = 4, dpi = 300)
      ggplot2::ggsave(reverse_path, reverse_plot, width = 8, height = 4, dpi = 300)
      tibble::tibble(
        sample_id = sample_id,
        forward_plot = forward_path,
        reverse_plot = reverse_path
      )
    }
  )
}
