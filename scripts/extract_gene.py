#!/usr/bin/env python
"""Extract sequences for a target gene from downloaded NCBI genomes.

Looks at both CDS FASTA (`*_cds_from_genomic.fna`) and rRNA FASTA
(`*_rna_from_genomic.fna`) files in the genomes directory, since rRNA
genes (e.g. 16S) are not in the CDS file.

Matching strategy:
  - Look at [gene=...] field (standard for CDS like tuf, rpoB)
  - Look at [product=...] field (standard for rRNA like "16S ribosomal RNA")
  - A genome contributes ONE sequence (the first match found), to avoid
    inflating the alignment with paralogs.

If no sequences are extracted, the script writes an empty FASTA and exits 0
with a warning. This lets the pipeline skip a gene that doesn't apply to
the chosen genus rather than crashing the whole run.
"""
import os
import re
import glob
import sys

from Bio import SeqIO


# Snakemake injects the `snakemake` object when run as a `script:`
genomes_dir = snakemake.input.genomes
output_path = snakemake.output[0]
aliases     = snakemake.params.aliases
gene_name   = snakemake.wildcards.gene

# normalize aliases to lowercase for case-insensitive matching
aliases = [a.lower() for a in aliases]

gene_re    = re.compile(r"\[gene=([^\]]+)\]",    re.IGNORECASE)
product_re = re.compile(r"\[product=([^\]]+)\]", re.IGNORECASE)


def header_matches(header: str) -> bool:
    """Return True if any alias matches the gene or product tag."""
    tags = []
    g = gene_re.search(header)
    if g:
        tags.append(g.group(1).lower())
    p = product_re.search(header)
    if p:
        tags.append(p.group(1).lower())
    if not tags:
        return False
    # match if any alias is a substring of any tag
    return any(any(a in tag for a in aliases) for tag in tags)


def candidate_files(genomes_dir: str):
    """Yield all CDS and rRNA FASTA files in the genomes directory."""
    patterns = [
        os.path.join(genomes_dir, "*_cds_from_genomic.fna"),
        os.path.join(genomes_dir, "*_rna_from_genomic.fna"),
    ]
    for pat in patterns:
        for f in glob.glob(pat):
            yield f


def genome_accession(filepath: str) -> str:
    """Strip suffix to get the genome accession (e.g. GCF_000009045.1)."""
    base = os.path.basename(filepath)
    for suffix in ("_cds_from_genomic.fna", "_rna_from_genomic.fna"):
        if base.endswith(suffix):
            return base[:-len(suffix)]
    return base


def main():
    written = 0
    skipped = 0

    # group files by genome accession (one CDS + one rRNA per genome)
    by_genome: dict[str, list[str]] = {}
    for f in candidate_files(genomes_dir):
        by_genome.setdefault(genome_accession(f), []).append(f)

    if not by_genome:
        sys.exit(f"[extract_gene/{gene_name}] No CDS or rRNA FASTA files "
                 f"found in {genomes_dir}.")

    with open(output_path, "w") as out:
        for accession, files in sorted(by_genome.items()):
            found = False
            for fna in files:
                if found:
                    break
                for rec in SeqIO.parse(fna, "fasta"):
                    if header_matches(rec.description):
                        rec.id = f"{accession}|{rec.id}"
                        rec.description = ""
                        SeqIO.write(rec, out, "fasta")
                        written += 1
                        found = True
                        break
            if not found:
                skipped += 1
                print(f"[extract_gene/{gene_name}] not found in {accession}",
                      file=sys.stderr)

    print(f"[extract_gene/{gene_name}] extracted {written} sequences, "
          f"{skipped} genomes had no match", file=sys.stderr)

    if written == 0:
        print(f"[extract_gene/{gene_name}] WARNING: no sequences extracted. "
              f"Check that '{gene_name}' (and its aliases) are valid for "
              f"this genus, or that the annotation uses these names.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
