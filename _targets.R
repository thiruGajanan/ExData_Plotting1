library(targets)

source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(source_files, source))

tar_option_set(
  packages = c(
    "readr", "dplyr", "tibble", "tidyr", "purrr", "stringr", "ggplot2",
    "ShortRead", "dada2", "decontam", "Biostrings", "DECIPHER", "phangorn",
    "phyloseq", "vegan", "ANCOMBC", "Maaslin2", "microbiome", "broom",
    "rlang", "biomformat"
  ),
  format = "rds"
)

config <- list(
  metadata_path = "data/metadata/sample_metadata.tsv",
  fastq_dir = "data/raw",
  reference_fasta = "references/silva_nr99_v138_train_set.fa.gz",
  species_fasta = "references/silva_species_assignment_v138.fa.gz",
  picrust2_threads = 4
)

list(
  tar_target(metadata_path, config$metadata_path, format = "file"),
  tar_target(metadata, load_metadata(metadata_path)),
  tar_target(fastq_dir, config$fastq_dir, format = "directory"),
  tar_target(fastq_table, find_fastqs(fastq_dir)),
  tar_target(quality_profiles, generate_quality_profiles(fastq_table), format = "rds"),
  tar_target(filtered_reads, filter_and_trim_reads(fastq_table)),
  tar_target(error_models, learn_error_models(filtered_reads)),
  tar_target(asv_results, infer_asvs(filtered_reads, error_models)),
  tar_target(seqtab_nochim, remove_chimeras(asv_results$sequence_table)),
  tar_target(contaminant_results, remove_contaminants(seqtab_nochim, metadata)),
  tar_target(taxonomy, assign_taxonomy(contaminant_results$seqtab, config$reference_fasta, config$species_fasta)),
  tar_target(phy_tree, build_phylogenetic_tree(contaminant_results$seqtab)),
  tar_target(physeq, make_phyloseq(contaminant_results$seqtab, taxonomy, metadata, phy_tree)),
  tar_target(alpha_diversity, compute_alpha_diversity(physeq)),
  tar_target(beta_diversity, compute_beta_diversity(physeq)),
  tar_target(ordinations, run_ordinations(physeq)),
  tar_target(permanova, run_permanova(beta_diversity$distance, metadata, stats::as.formula("distance ~ condition + batch"))),
  tar_target(ancombc_results, run_ancombc(physeq, formula = "condition + batch")),
  tar_target(maaslin2_results, run_maaslin2(physeq, metadata, fixed_effects = c("condition", "batch"))),
  tar_target(core_microbiome, summarize_core_microbiome(physeq)),
  tar_target(function_predictions, predict_functional_profiles(physeq, method = "picrust2", threads = config$picrust2_threads)),
  tar_target(pathway_analysis, analyze_pathway_differences(function_predictions$pathway_abundance, metadata, group = condition))
)
