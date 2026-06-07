#!/usr/bin/env python3
"""db_gene_counts.py - how many marker genes each species' reference genome has in
the protal DB.

protal's gene universe is a fixed set of universal marker genes (gene2geneid.tsv,
168 in r226). A species' representative genome contains only a subset of them.
This count is the true denominator for "how many genes could possibly appear" --
the strain meta.tsv only lists genes that got >=1 read hit, so comparing the two
reveals how many markers were lost to abundance/coverage *before* M3.

Maps each strain species -> its genome internal id (genome2tiid.tsv) and counts
that genome's rows in unique_kmers.tsv. Writes <out> as:
    species <tab> genome_id <tab> db_markers

Usage:
  db_gene_counts.py --db <DB_DIR> --strains <strains_dir> --out <db_gene_counts.tsv>
"""
import argparse
import glob
import os
import sys


def norm(name):
    """Normalise a species label for matching ('s__Foo bar' <-> 's__Foo_bar')."""
    return name.replace(" ", "_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="protal DB directory")
    ap.add_argument("--strains", required=True, help="strain output dir (for species list)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    species = [norm(os.path.basename(m)[:-len(".meta.tsv")])
               for m in glob.glob(os.path.join(args.strains, "*.meta.tsv"))]
    want = set(species)

    # species -> set of genome internal ids (last taxonomy field is s__...)
    g2t = os.path.join(args.db, "genome2tiid.tsv")
    sp_to_ids = {}
    with open(g2t, errors="ignore") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            gid = f[1]
            sp = norm(f[3].split(";")[-1].strip())
            if sp in want:
                sp_to_ids.setdefault(sp, set()).add(gid)
    id_to_sp = {gid: sp for sp, ids in sp_to_ids.items() for gid in ids}

    # one streaming pass over the big unique_kmers.tsv, counting rows per wanted id
    uk = os.path.join(args.db, "unique_kmers.tsv")
    per_id = {gid: 0 for gid in id_to_sp}
    with open(uk, errors="ignore") as fh:
        for line in fh:
            i = line.find("\t")
            if i < 0:
                continue
            gid = line[:i]
            if gid in per_id:
                per_id[gid] += 1

    # per species: take the genome with the most markers (the representative)
    with open(args.out, "w") as out:
        out.write("species\tgenome_id\tdb_markers\n")
        for sp in sorted(want):
            ids = sp_to_ids.get(sp)
            if not ids:
                out.write(f"{sp}\t\t\n")
                continue
            best = max(ids, key=lambda g: per_id.get(g, 0))
            out.write(f"{sp}\t{best}\t{per_id.get(best, 0)}\n")
    sys.stderr.write(f"Wrote {args.out} ({len(want)} species)\n")


if __name__ == "__main__":
    main()
