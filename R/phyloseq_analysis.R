#' Create a phyloseq object from ASV data
#'
#' @param seqtab_clean Sequence table with contaminants removed.
#' @param taxonomy Taxonomy table returned by `assign_taxonomy`.
#' @param metadata Tibble of sample metadata.
#' @param phy_tree A phylogenetic tree object.
#' @return A `phyloseq` object bundling counts, taxonomy, metadata, and tree.
make_phyloseq <- function(seqtab_clean, taxonomy, metadata, phy_tree) {
  otu <- phyloseq::otu_table(seqtab_clean, taxa_are_rows = FALSE)
  tax <- phyloseq::tax_table(as.matrix(taxonomy))
  rownames(metadata) <- metadata$sample_id
  samp <- phyloseq::sample_data(metadata)
  phyloseq::phyloseq(otu, tax, samp, phy_tree = phy_tree)
}

#' Compute alpha-diversity metrics
#'
#' @param physeq A `phyloseq` object.
#' @param measures Vector of alpha-diversity measures to calculate.
#' @return Tibble containing alpha-diversity values per sample.
compute_alpha_diversity <- function(physeq, measures = c("Observed", "Shannon", "Simpson")) {
  richness <- phyloseq::estimate_richness(physeq, measures = measures)
  tibble::as_tibble(richness, rownames = "sample_id")
}

#' Compute beta-diversity distance matrices
#'
#' @param physeq A `phyloseq` object.
#' @param distance Metric passed to `phyloseq::distance`.
#' @return A list containing the distance object and a tidy representation.
compute_beta_diversity <- function(physeq, distance = "bray") {
  dist_obj <- phyloseq::distance(physeq, method = distance)
  dist_tbl <- as.matrix(dist_obj) |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tidyr::pivot_longer(-sample_id, names_to = "comparison", values_to = "distance")
  list(distance = dist_obj, distance_tidy = dist_tbl)
}

#' Ordinate a phyloseq object using multiple methods
#'
#' @param physeq A `phyloseq` object.
#' @param methods Ordination techniques to evaluate.
#' @param distance Distance metric shared across ordinations.
#' @return A list of ordination objects and tidy coordinates.
run_ordinations <- function(physeq,
                            methods = c("PCoA", "NMDS"),
                            distance = "bray") {
  purrr::map(methods, function(method) {
    ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
    coords <- phyloseq::plot_ordination(physeq, ord, justDF = TRUE)
    list(method = method, ordination = ord, coordinates = coords)
  }) |>
    stats::setNames(methods)
}

#' Perform PERMANOVA for group comparisons
#'
#' @param distance_obj A distance object from `compute_beta_diversity`.
#' @param metadata Sample metadata tibble.
#' @param formula Model formula passed to `vegan::adonis2`.
#' @param permutations Number of permutations.
#' @return Tidy tibble of PERMANOVA results.
run_permanova <- function(distance_obj,
                          metadata,
                          formula,
                          permutations = 999) {
  sample_ids <- attr(distance_obj, "Labels")
  if (is.null(sample_ids)) {
    sample_ids <- rownames(as.matrix(distance_obj))
  }
  metadata_aligned <- metadata |>
    dplyr::filter(sample_id %in% sample_ids) |>
    dplyr::arrange(match(sample_id, sample_ids))
  if (nrow(metadata_aligned) == 0) {
    stop("Metadata does not contain the samples present in the distance object.")
  }
  rownames(metadata_aligned) <- metadata_aligned$sample_id
  distance <- distance_obj
  permanova <- vegan::adonis2(formula, data = metadata_aligned, permutations = permutations, by = "margin")
  permanova |>
    broom::tidy() |>
    tibble::as_tibble()
}

#' Summarise ordination coordinates with metadata for plotting
#'
#' @param ordinations Output from `run_ordinations`.
#' @param metadata Sample metadata tibble.
#' @return Tibble with ordination axes joined to metadata.
combine_ord_coord <- function(ordinations, metadata) {
  purrr::imap_dfr(
    ordinations,
    function(obj, method) {
      tibble::as_tibble(obj$coordinates) |>
        dplyr::rename(sample_id = SampleID) |>
        dplyr::left_join(metadata, by = "sample_id") |>
        dplyr::mutate(method = method)
    }
  )
}
