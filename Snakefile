# ============================================================================
# genus-primer-pipeline
#
# For a given bacterial genus and one or more housekeeping genes:
#   download genomes -> extract gene -> dereplicate -> align ->
#   design primers -> in silico PCR validation -> HTML report
#
# Usage:
#   snakemake --use-conda --cores 4
# ============================================================================

configfile: "config.yaml"

wildcard_constraints:
     gene = r"[^./]+"

GENUS = config["genus"]
GENES = config["genes"]

# Per-genus result tree keeps multiple genus runs cleanly separated
RESULTS   = f"results/{GENUS}"
GENOMES   = f"{RESULTS}/genomes"
EXTRACTED = f"{RESULTS}/extracted"
ALIGNED   = f"{RESULTS}/aligned"
PRIMERS   = f"{RESULTS}/primers"


# ---------------------------------------------------------------------------
# Top-level target
# ---------------------------------------------------------------------------
rule all:
    input:
        expand(f"{PRIMERS}/{{gene}}_report.html", gene=GENES),
        # only build comparison report when multiple genes are configured
        f"{RESULTS}/comparison.html" if len(GENES) > 1 else []


# ---------------------------------------------------------------------------
# 1. Resolve genus name -> NCBI taxid
# ---------------------------------------------------------------------------
rule resolve_taxid:
    output:
        f"{RESULTS}/taxid.txt"
    params:
        genus=GENUS,
        email=config["ncbi_email"]
    conda:
        "envs/entrez.yaml"
    shell:
        "python scripts/resolve_taxid.py "
        "--genus {params.genus} --email {params.email} --output {output}"


# ---------------------------------------------------------------------------
# 2. Download all complete genomes for the genus
# ---------------------------------------------------------------------------
checkpoint download_genomes:
    input:
        f"{RESULTS}/taxid.txt"
    output:
        directory(GENOMES)
    params:
        domain=config["domain"]
    conda:
        "envs/ncbi_download.yaml"
    shell:
        r"""
        TAXID=$(cat {input})
        mkdir -p {output}
        ncbi-genome-download {params.domain} \
            -F 'cds-fasta,fasta,rna-fasta' \
            -l 'complete' \
            --species-taxids "$TAXID" \
            --flat-output \
            -o {output}
        find {output} -name "*.gz" -exec gunzip {{}} \;
        """


# ---------------------------------------------------------------------------
# 3. Extract target gene from each genome
# ---------------------------------------------------------------------------
rule extract_gene:
    input:
        genomes=GENOMES
    output:
        f"{EXTRACTED}/{{gene}}.fasta"
    params:
        aliases=lambda wc: config.get("gene_aliases", {}).get(wc.gene, [wc.gene])
    conda:
        "envs/seqtools.yaml"
    script:
        "scripts/extract_gene.py"


# ---------------------------------------------------------------------------
# 4. Dereplicate identical sequences (keeps centroids)
# ---------------------------------------------------------------------------
rule cluster:
    input:
        f"{EXTRACTED}/{{gene}}.fasta"
    output:
        f"{EXTRACTED}/{{gene}}.centroids.fasta"
    conda:
        "envs/seqtools.yaml"
    shell:
        "vsearch --cluster_fast {input} --strand both "
        "--id 0.97 --centroids {output} --quiet"


# ---------------------------------------------------------------------------
# 5. Multiple sequence alignment
# ---------------------------------------------------------------------------
rule align:
    input:
        f"{EXTRACTED}/{{gene}}.centroids.fasta"
    output:
        f"{ALIGNED}/{{gene}}.aln"
    conda:
        "envs/seqtools.yaml"
    shell:
        "muscle -align {input} -output {output} 2> /dev/null"


# ---------------------------------------------------------------------------
# 6. Scan alignment for primer candidates (R)
# ---------------------------------------------------------------------------
rule design_primers:
    input:
        alignment=f"{ALIGNED}/{{gene}}.aln"
    output:
        primers=f"{PRIMERS}/{{gene}}.primers.tsv",
        diversity_plot=f"{PRIMERS}/{{gene}}.diversity.png"
    params:
        cfg=config["primer"]
    conda:
        "envs/r.yaml"
    script:
        "scripts/design_primers.R"


# ---------------------------------------------------------------------------
# 7. In silico PCR validation against all genomes
# ---------------------------------------------------------------------------
rule in_silico_pcr:
    input:
        primers=f"{PRIMERS}/{{gene}}.primers.tsv",
        genomes=GENOMES
    output:
        f"{PRIMERS}/{{gene}}.pcr_results.tsv"
    conda:
        "envs/seqtools.yaml"
    script:
        "scripts/in_silico_pcr.py"


# ---------------------------------------------------------------------------
# 8a. Per-gene HTML report
# ---------------------------------------------------------------------------
rule gene_report:
    input:
        primers=f"{PRIMERS}/{{gene}}.primers.tsv",
        pcr=f"{PRIMERS}/{{gene}}.pcr_results.tsv",
        plot=f"{PRIMERS}/{{gene}}.diversity.png"
    output:
        f"{PRIMERS}/{{gene}}_report.html"
    params:
        genus=GENUS
    conda:
        "envs/r.yaml"
    script:
        "scripts/gene_report.Rmd"


# ---------------------------------------------------------------------------
# 8b. Cross-gene comparison report (only when len(GENES) > 1)
# ---------------------------------------------------------------------------
rule comparison_report:
    input:
        primers=expand(f"{PRIMERS}/{{gene}}.primers.tsv", gene=GENES),
        pcr=expand(f"{PRIMERS}/{{gene}}.pcr_results.tsv", gene=GENES)
    output:
        f"{RESULTS}/comparison.html"
    params:
        genus=GENUS,
        genes=GENES
    conda:
        "envs/r.yaml"
    script:
        "scripts/comparison_report.Rmd"
