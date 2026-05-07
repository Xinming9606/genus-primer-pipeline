#!/usr/bin/env python
"""Run in silico PCR using `seqkit amplicon` for each candidate primer pair.

For each pair listed in the primers TSV, attempts amplification against every
genome FASTA (`*_genomic.fna`, excluding the CDS / rRNA derived files). Reports
per-pair retrieval rate (fraction of genomes with >=1 hit) and amplicon length
statistics. Tolerates up to 3 mismatches per primer (configurable below).
"""
import csv
import glob
import os
import statistics
import subprocess
import sys
import tempfile


# Snakemake injection
primers_tsv  = snakemake.input.primers
genomes_dir  = snakemake.input.genomes
output_path  = snakemake.output[0]

MAX_MISMATCHES = 3   # passed to seqkit amplicon -m


def find_full_genomes(genomes_dir: str):
    """Return whole-genome FASTA files, excluding CDS/RNA derived files."""
    all_fna = glob.glob(os.path.join(genomes_dir, "*_genomic.fna"))
    return [
        f for f in all_fna
        if not f.endswith("_cds_from_genomic.fna")
        and not f.endswith("_rna_from_genomic.fna")
    ]


def parse_amplicons(fasta_text: str):
    """Yield amplicon lengths from a FASTA-formatted seqkit output string."""
    if not fasta_text.strip():
        return
    current_seq = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if current_seq:
                yield len("".join(current_seq))
                current_seq = []
        else:
            current_seq.append(line.strip())
    if current_seq:
        yield len("".join(current_seq))


def main():
    if not os.path.exists(primers_tsv) or os.path.getsize(primers_tsv) == 0:
        print(f"[in_silico_pcr] primers TSV is empty; writing empty output",
              file=sys.stderr)
        with open(output_path, "w") as out:
            out.write("primer_id\tfwd\trev\tn_hits\tn_genomes\t"
                      "retrieval_rate\tmean_amplicon_len\tsd_amplicon_len\n")
        return

    genomes = find_full_genomes(genomes_dir)
    n_genomes = len(genomes)
    if n_genomes == 0:
        sys.exit(f"[in_silico_pcr] no whole-genome FASTA files found in "
                 f"{genomes_dir}")

    print(f"[in_silico_pcr] {n_genomes} genomes available", file=sys.stderr)

    with open(primers_tsv) as f, open(output_path, "w") as out:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "primer_id", "fwd", "rev",
            "n_hits", "n_genomes", "retrieval_rate",
            "mean_amplicon_len", "sd_amplicon_len"
        ])

        for row in reader:
            pid = row["primer_id"]
            fwd = row["fwd"]
            rev = row["rev"]

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False
            ) as pf:
                pf.write(f"{pid}\t{fwd}\t{rev}\n")
                primer_file = pf.name

            try:
                hits = 0
                lengths = []
                for genome in genomes:
                    try:
                        res = subprocess.run(
                            ["seqkit", "amplicon",
                             "-p", primer_file,
                             "-m", str(MAX_MISMATCHES),
                             genome],
                            capture_output=True, text=True, check=True,
                        )
                    except subprocess.CalledProcessError as e:
                        print(f"[in_silico_pcr] seqkit failed on "
                              f"{os.path.basename(genome)} for {pid}: "
                              f"{e.stderr.strip()}", file=sys.stderr)
                        continue

                    amplicon_lengths = list(parse_amplicons(res.stdout))
                    if amplicon_lengths:
                        hits += 1
                        lengths.extend(amplicon_lengths)
            finally:
                os.unlink(primer_file)

            rate = hits / n_genomes if n_genomes else 0
            mean_len = statistics.mean(lengths) if lengths else 0
            sd_len = statistics.stdev(lengths) if len(lengths) > 1 else 0

            writer.writerow([
                pid, fwd, rev,
                hits, n_genomes, f"{rate:.3f}",
                f"{mean_len:.1f}", f"{sd_len:.1f}"
            ])
            print(f"[in_silico_pcr] {pid}: {hits}/{n_genomes} genomes "
                  f"({rate:.1%})", file=sys.stderr)


if __name__ == "__main__":
    main()
