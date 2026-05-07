# genus-primer-pipeline

A Snakemake pipeline for designing amplicon-sequencing primers for any
bacterial genus and housekeeping gene of interest.

## What it does

Given a genus name (e.g. `Bacillus`) and one or more housekeeping genes
(e.g. `tuf`, `rpoB`), the pipeline:

1. Resolves the genus name to an NCBI taxid
2. Downloads all complete genomes for that genus
3. Extracts the target gene from each genome (handles both CDS and rRNA)
4. Aligns the sequences (MUSCLE) after dereplication (vsearch)
5. Scans the alignment for conserved primer sites flanking variable regions
6. Validates candidate primer pairs by in silico PCR
7. Outputs an HTML report ranking the top primer pairs

When multiple genes are specified, an additional comparison report is
generated.

## Installation

The pipeline needs **conda** (or **mamba** / **micromamba**) and an active
internet connection on first run. All other dependencies are installed
automatically into per-rule environments by Snakemake.

```bash
# 1. Clone
git clone git@github.com:Xinming9606/genus-primer-pipeline.git
cd genus-primer-pipeline

# 2. Create the runner environment (only needs snakemake itself)
conda env create -f envs/runner.yaml
conda activate primer-pipeline
```

## Usage

```bash
# 1. Edit config.yaml: set genus, genes, ncbi_email
# 2. Dry-run first to verify the DAG
snakemake --dry-run

# 3. Real run (--use-conda creates per-rule envs the first time)
snakemake --use-conda --cores 4
```

Outputs land under `results/<genus>/`. The main deliverable is
`results/<genus>/primers/<gene>_report.html` (per gene), plus
`results/<genus>/comparison.html` if more than one gene is configured.

## Repository layout
