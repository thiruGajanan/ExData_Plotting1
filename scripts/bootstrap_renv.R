# Run this script to (re)initialise the renv environment programmatically.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
if (file.exists("renv.lock")) {
  renv::activate()
} else {
  renv::init(bare = TRUE)
}
packages <- c(
  "targets", "dada2", "ShortRead", "decontam", "Biostrings", "DECIPHER", "phyloseq",
  "phangorn", "vegan", "ANCOMBC", "Maaslin2", "microbiome", "biomformat",
  "ggplot2", "dplyr", "tidyr", "tibble", "readr", "purrr", "stringr", "broom", "rlang"
)
renv::install(packages)
renv::snapshot(prompt = FALSE)
