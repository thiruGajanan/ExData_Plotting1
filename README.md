# Amplicon Sequencing Workflow

This repository implements an end-to-end 16S rRNA gene amplicon workflow that imports raw FASTQ files, infers amplicon sequence variants (ASVs) with **dada2**, removes contaminants with **decontam**, constructs a **phyloseq** object, performs diversity and ordination analyses, tests for differential abundance, predicts functional potential, and documents the results in R Markdown. The pipeline is orchestrated with **targets** and dependencies are tracked with **renv** for full reproducibility.

## Repository structure

```
analysis/                      R Markdown report that summarises pipeline outputs
R/                             Modular R functions used by the pipeline
_targets.R                     targets pipeline definition
renv.lock                      Locked package versions for reproducibility
data/
  raw/                         Place raw paired-end FASTQ files here
  metadata/                    Sample metadata table (TSV)
outputs/                       Derived figures and tables
references/                    Reference training sets (e.g., SILVA/GTDB)
scripts/                       Helper scripts and automation entry points
```

Existing plotting scripts from the original coursework submission (`plot1.R` – `plot4.R`) are retained for reference but are not part of the new microbiome workflow.

## Getting started

1. **Install dependencies with renv**
   ```r
   install.packages("renv")
   renv::restore()
   ```

2. **Prepare inputs**
   - Copy paired-end FASTQ files into `data/raw/` with names such as `SampleA_R1_001.fastq.gz` / `SampleA_R2_001.fastq.gz`.
   - Update `data/metadata/sample_metadata.tsv` with sample-specific information (include an `is_control` column to flag negative controls).
   - Download SILVA or GTDB reference training sets into `references/` and adjust the paths in `_targets.R` if required.

3. **Run the pipeline**
   ```r
   renv::activate()
   targets::tar_make()
   ```

4. **Render the report**
   ```r
   rmarkdown::render("analysis/microbiome_workflow.Rmd")
   ```

## Workflow summary

1. **Import FASTQ files and metadata** using `ShortRead` and `dada2`, producing read quality diagnostics saved under `outputs/quality_profiles/`.
2. **Trim, filter, and infer ASVs** with `dada2`, learn error models, merge paired reads, remove chimeras, and remove contaminants (`decontam`).
3. **Assign taxonomy** (SILVA/GTDB) and build a comprehensive `phyloseq` object linking ASV counts, taxonomy, sample metadata, and an inferred phylogenetic tree.
4. **Characterise alpha/beta diversity**, generate ordinations, and test for group-level differences via PERMANOVA (`vegan`).
5. **Perform differential abundance testing** with ANCOM-BC and MaAsLin2, complemented by core microbiome summaries (`microbiome`).
6. **Predict functional profiles** via PICRUSt2 (or optionally Tax4Fun2) and highlight pathway shifts between experimental groups.
7. **Document the analysis** in `analysis/microbiome_workflow.Rmd`, leveraging `targets` for reproducible execution and `renv` for environment capture.

## Key scripts

- `_targets.R` — defines the computational pipeline and dependencies.
- `analysis/microbiome_workflow.Rmd` — narrative report consuming pipeline outputs.
- `R/*.R` — modular functions for data import, ASV processing, phyloseq analyses, differential abundance, and functional prediction.

## Notes

- The repository includes template metadata and directory placeholders but does not ship with sequencing data or reference databases.
- Update the `config` list in `_targets.R` to reflect project-specific parameters (e.g., filtering thresholds, PICRUSt2 threads, reference file locations).
- PICRUSt2 must be installed separately and accessible on the system `PATH` (see [https://github.com/picrust/picrust2](https://github.com/picrust/picrust2)). If unavailable, the functional prediction step returns informative warnings without halting the pipeline.
