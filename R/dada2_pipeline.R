#' Filter and trim paired-end reads
#'
#' @param fastq_table Tibble produced by `find_fastqs`.
#' @param filt_dir Output directory for filtered reads.
#' @param trunc_len Vector of length two with truncation lengths for forward and reverse reads.
#' @param max_ee Vector of expected error thresholds.
#' @param trunc_q Quality score threshold for truncation.
#' @param max_n Maximum number of ambiguous bases allowed.
#' @param multithread Use multithreading where available.
#' @return Tibble describing filtered FASTQ files and read counts.
filter_and_trim_reads <- function(fastq_table,
                                  filt_dir = file.path("data", "filtered"),
                                  trunc_len = c(240, 200),
                                  max_ee = c(2, 2),
                                  trunc_q = 2,
                                  max_n = 0,
                                  min_len = 50,
                                  rm_phix = TRUE,
                                  compress = TRUE,
                                  multithread = TRUE) {
  if (nrow(fastq_table) == 0) {
    stop("No FASTQ files supplied to filter.")
  }
  dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
  filt_forward <- file.path(filt_dir, paste0(fastq_table$sample_id, "_F_filt.fastq.gz"))
  filt_reverse <- file.path(filt_dir, paste0(fastq_table$sample_id, "_R_filt.fastq.gz"))
  out <- dada2::filterAndTrim(
    fastq_table$forward,
    filt_forward,
    fastq_table$reverse,
    filt_reverse,
    truncLen = trunc_len,
    maxEE = max_ee,
    truncQ = trunc_q,
    maxN = max_n,
    rm.phix = rm_phix,
    compress = compress,
    multithread = multithread,
    minLen = min_len
  )
  tibble::tibble(
    sample_id = fastq_table$sample_id,
    filt_forward = filt_forward,
    filt_reverse = filt_reverse,
    reads_in = out[, 1],
    reads_out = out[, 2]
  )
}

#' Learn DADA2 error models for forward and reverse reads
#'
#' @param filtered_reads Tibble returned by `filter_and_trim_reads`.
#' @param multithread Use multithreading where available.
#' @return A list with forward and reverse error models.
learn_error_models <- function(filtered_reads, multithread = TRUE) {
  errF <- dada2::learnErrors(filtered_reads$filt_forward, multithread = multithread)
  errR <- dada2::learnErrors(filtered_reads$filt_reverse, multithread = multithread)
  list(errF = errF, errR = errR)
}

#' Infer ASVs using the DADA2 core algorithm
#'
#' @param filtered_reads Tibble returned by `filter_and_trim_reads`.
#' @param error_models Output from `learn_error_models`.
#' @param pooling Pooling strategy passed to `dada2::dada`.
#' @return A list containing the sequence table, mergers, and dereplication objects.
infer_asvs <- function(filtered_reads,
                       error_models,
                       pooling = "pseudo",
                       multithread = TRUE) {
  derep_forward <- purrr::map(filtered_reads$filt_forward, dada2::derepFastq)
  derep_reverse <- purrr::map(filtered_reads$filt_reverse, dada2::derepFastq)
  names(derep_forward) <- filtered_reads$sample_id
  names(derep_reverse) <- filtered_reads$sample_id
  dada_forward <- purrr::map(derep_forward, dada2::dada, err = error_models$errF, pool = pooling, multithread = multithread)
  dada_reverse <- purrr::map(derep_reverse, dada2::dada, err = error_models$errR, pool = pooling, multithread = multithread)
  mergers <- purrr::pmap(
    list(dada_forward, derep_forward, dada_reverse, derep_reverse),
    dada2::mergePairs,
    verbose = TRUE
  )
  sequence_table <- dada2::makeSequenceTable(mergers)
  list(
    sequence_table = sequence_table,
    mergers = mergers,
    derep_forward = derep_forward,
    derep_reverse = derep_reverse,
    dada_forward = dada_forward,
    dada_reverse = dada_reverse
  )
}

#' Remove chimeric sequences from a sequence table
#'
#' @param sequence_table Sequence table returned by `infer_asvs`.
#' @return A chimera-filtered sequence table.
remove_chimeras <- function(sequence_table) {
  dada2::removeBimeraDenovo(sequence_table, method = "consensus", multithread = TRUE, verbose = TRUE)
}

#' Identify and remove contaminant ASVs
#'
#' @param seqtab_nochim Chimera-filtered sequence table.
#' @param metadata Sample metadata tibble.
#' @param control_column Column indicating negative controls (logical or 0/1).
#' @param method Contaminant detection method passed to `decontam::isContaminant`.
#' @return A list containing the cleaned sequence table and contaminant statistics.
remove_contaminants <- function(seqtab_nochim,
                                metadata,
                                control_column = "is_control",
                                method = "prevalence") {
  if (!control_column %in% names(metadata)) {
    stop("Control column ", control_column, " not found in metadata.")
  }
  common_samples <- intersect(rownames(seqtab_nochim), metadata$sample_id)
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between sequence table and metadata.")
  }
  seqtab_aligned <- seqtab_nochim[common_samples, , drop = FALSE]
  metadata_aligned <- metadata |>
    dplyr::filter(sample_id %in% common_samples) |>
    dplyr::arrange(match(sample_id, common_samples))
  neg <- metadata_aligned[[control_column]]
  if (!is.logical(neg)) {
    neg <- as.logical(neg)
  }
  if (any(is.na(neg))) {
    stop("Negative control column contains missing values.")
  }
  rownames(seqtab_aligned) <- metadata_aligned$sample_id
  contam <- decontam::isContaminant(seqtab_aligned, neg = neg, method = method)
  seqtab_clean <- seqtab_aligned[, !contam$contaminant, drop = FALSE]
  list(
    seqtab = seqtab_clean,
    contaminants = tibble::as_tibble(contam, rownames = "asv"),
    removed_asvs = tibble::as_tibble(contam, rownames = "asv") |>
      dplyr::filter(contaminant)
  )
}

#' Assign taxonomy to ASVs using reference training sets
#'
#' @param seqtab_clean Contaminant-filtered sequence table.
#' @param reference_fasta Path to the training set FASTA (e.g., SILVA or GTDB).
#' @param species_fasta Optional species-level training set.
#' @return Taxonomy matrix compatible with `phyloseq::tax_table`.
assign_taxonomy <- function(seqtab_clean,
                            reference_fasta,
                            species_fasta = NULL,
                            multithread = TRUE) {
  if (!file.exists(reference_fasta)) {
    stop("Reference FASTA not found: ", reference_fasta)
  }
  sequences <- colnames(seqtab_clean)
  taxa <- dada2::assignTaxonomy(sequences, reference_fasta, multithread = multithread)
  if (!is.null(species_fasta)) {
    taxa <- dada2::addSpecies(taxa, species_fasta, allowMultiple = TRUE)
  }
  rownames(taxa) <- sequences
  taxa
}

#' Construct a phylogenetic tree from ASV sequences
#'
#' @param seqtab_clean Contaminant-filtered sequence table.
#' @return A rooted phylogenetic tree.
build_phylogenetic_tree <- function(seqtab_clean) {
  sequences <- Biostrings::DNAStringSet(colnames(seqtab_clean))
  alignment <- DECIPHER::AlignSeqs(sequences, processors = NULL)
  phangorn::phyDat(alignment, type = "DNA") |>
    {
      dm <- phangorn::dist.ml(., model = "JC")
      treeNJ <- phangorn::NJ(dm)
      fit <- phangorn::pml(treeNJ, data = .)
      fit_opt <- phangorn::optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")
      phangorn::midpoint(fit_opt$tree)
    }
}
