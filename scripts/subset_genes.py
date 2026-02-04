#!/usr/bin/env python3
"""
Subset a reference FASTA by a newline-delimited gene list.

Header format supported:
  >taxid_geneid
  >taxid_geneid_COMMENT    # comment allowed after the id

Usage:
  subset_genes.py --ref REF.fasta --build_gene_subset genes.txt --out subset.fasta

Options:
  --build_gene_subset  newline-delimited gene list (taxid_geneid)

This script purposely avoids third-party deps and streams the FASTA.
"""

from __future__ import annotations
import argparse
import gzip
import sys
from typing import Set, TextIO, Tuple


def open_maybe_gz(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_gene_id_from_header(header_token: str) -> Tuple[str, str]:
    """Return (taxid_geneid, gene_only) from a header token.

    Examples:
      '2463_11' -> ('2463_11', '11')
      '2463_11_comment' -> ('2463_11', '11')
      '11' -> ('11', '11')
    """
    parts = header_token.split("_")
    if len(parts) >= 2:
        return (f"{parts[0]}_{parts[1]}", parts[1])
    return (header_token, header_token)


def parse_gene_list(path: str) -> Tuple[Set[str], Set[str]]:
    """Return (full_ids, gene_only_ids).

    full_ids contain entries like 'taxid_geneid'.
    gene_only_ids contain entries like '11' to match any taxid.
    """
    full: Set[str] = set()
    gene_only: Set[str] = set()
    with open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            token = line.split()[0]
            # token like 'taxid_geneid' or 'taxid' 'geneid'
            tokens_under = token.split("_")
            if len(tokens_under) >= 2:
                full.add(f"{tokens_under[0]}_{tokens_under[1]}")
            else:
                sp = line.split()
                if len(sp) >= 2:
                    full.add(f"{sp[0]}_{sp[1]}")
                else:
                    # treat single token as gene-only id
                    gene_only.add(token)
    return full, gene_only


def subset_fasta(ref: str, full_genes: Set[str], gene_only: Set[str], out: str):
    total = 0
    written = 0

    with open_maybe_gz(ref, "rt") as infh, open_maybe_gz(out, "wt") as outfh:
        header = None
        seq_lines = []
        for line in infh:
            if line.startswith(">"):
                if header is not None:
                    total += 1
                    token = header[1:].strip().split()[0]
                    taxid_geneid, gene_only_id = normalize_gene_id_from_header(token)
                    selected = False
                    if taxid_geneid in full_genes:
                        selected = True
                    if gene_only_id in gene_only:
                        selected = True

                    if selected:
                        outfh.write(header)
                        outfh.writelines(seq_lines)
                        written += 1

                header = line
                seq_lines = []
            else:
                if header is None:
                    continue
                seq_lines.append(line)

        if header is not None:
            total += 1
            token = header[1:].strip().split()[0]
            taxid_geneid, gene_only_id = normalize_gene_id_from_header(token)
            selected = False
            if taxid_geneid in full_genes:
                selected = True
            if gene_only_id in gene_only:
                selected = True
            if selected:
                outfh.write(header)
                outfh.writelines(seq_lines)
                written += 1

    # Do not track or report missing IDs; caller can infer from counts if needed
    return total, written


def main(argv=None):
    p = argparse.ArgumentParser(description="Subset a FASTA by gene list (taxid_geneid)")
    p.add_argument("--ref", required=True, help="reference FASTA (can be .gz)")
    p.add_argument("--build_gene_subset", required=True, help="newline-delimited gene list (taxid_geneid or gene-only)")
    p.add_argument("--out", required=True, help="output FASTA (can be .gz)")
    args = p.parse_args(argv)

    full_genes, gene_only = parse_gene_list(args.build_gene_subset)
    if not full_genes and not gene_only:
        print("No genes parsed from gene list; aborting.", file=sys.stderr)
        return 2

    total, written = subset_fasta(args.ref, full_genes, gene_only, args.out)

    requested = len(full_genes) + len(gene_only)

    print(f"Total records scanned: {total}")
    print(f"Records written: {written}")
    print(f"Genes requested: {requested}")

    # We intentionally do not track or report missing IDs; script exits 0 on success
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
