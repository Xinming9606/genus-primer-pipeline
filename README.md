# genus-primer-pipeline

A Snakemake pipeline for designing amplicon-sequencing primers for any
bacterial genus and housekeeping gene of interest.

## What it does

Given a genus name (e.g. `Bacillus`) and one or more housekeeping genes
(e.g. `tuf`, `rpoB`), the pipeline:

1. Resolves the genus name to an NCBI taxid
2. Downloads all complete genomes for that genus
3. Extracts the target gene from each genome
4. Aligns the sequences (MUSCLE) after dereplication (vsearch)
5. Scans the alignment for conserved primer sites flanking variable regions
6. Validates candidate primer pairs by in silico PCR
7. Outputs an HTML report ranking the top primer pairs

When multiple genes are specified, an additional comparison report is
generated.

## Quick start

```bash
# 1. Edit config.yaml: set genus, genes, ncbi_email
# 2. Run
snakemake --use-conda --cores 4
```

## Status

Work in progress. See commit history for development trajectory.
