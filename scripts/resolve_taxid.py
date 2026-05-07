#!/usr/bin/env python
"""Resolve a genus name to an NCBI taxid via Entrez.

Writes the taxid as a single line to the output file. Fails loudly if the
genus is not found or is ambiguous, so the pipeline doesn't silently
download the wrong organisms.
"""
import argparse
import sys
from Bio import Entrez


def resolve(genus: str, email: str) -> str:
    Entrez.email = email
    handle = Entrez.esearch(db="taxonomy", term=f"{genus}[Scientific Name]")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        sys.exit(f"[resolve_taxid] No taxid found for genus '{genus}'. "
                 f"Check spelling and that it is a valid scientific name.")
    if len(ids) > 1:
        sys.exit(f"[resolve_taxid] Ambiguous match for '{genus}' "
                 f"(got {len(ids)} taxids: {ids}). Use a more specific name.")
    return ids[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genus", required=True)
    ap.add_argument("--email", required=True,
                    help="Email registered with NCBI Entrez (required by API).")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    taxid = resolve(args.genus, args.email)
    with open(args.output, "w") as f:
        f.write(taxid + "\n")
    print(f"[resolve_taxid] {args.genus} -> taxid {taxid}", file=sys.stderr)


if __name__ == "__main__":
    main()
