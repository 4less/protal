#!/usr/bin/env python3
"""qcmsa.py - QC and filter a strain MSA using protal .meta.tsv quality metrics.

This is the Python port of qcmsa.R (milestone M5). It is dependency-free for the
core filtering (Python >= 3.8 standard library only); matplotlib is imported
lazily and only when --plot is requested.

The filter has two passes:

  1. MRate2 (multi-allelicity) filter  -- identical algorithm to qcmsa.R.
     Iterative Tukey-IQR outlier detection on the *count* of multi-allelic
     (MRate2 > 0) cells per gene and per sample. Because MRate2 is sparse
     (mostly 0) a plain fence on the rates collapses to 0 and fires on any
     non-zero value, so we fence the non-zero counts instead and additionally
     require at least --min-bad bad peers before anything is removed.

  2. Site / sequence cleanup (M5 step 4b) -- optional second pass over the
     surviving columns: drop constant and near-constant (low-parsimony) sites,
     optionally mask individual cell outliers, and optionally re-apply the
     per-sequence horizontal-coverage floor after gene removal.

Inputs match protal's output contract:
  <msa>        FASTA (plain or .gz) -- protal's <species>.raw.msa.fna
  <partition>  RAxML-style partition -- "DNA, gene<ID> = <start>-<end>" (0-based inclusive)
  <meta.tsv>   protal per-sample x per-gene metrics, WITH a header row

Usage:
  qcmsa.py <msa> <partition> <meta.tsv> [options]
"""

import argparse
import gzip
import math
import os
import re
import sys
from collections import Counter, defaultdict

# Characters treated as "missing" (no informative base) at an MSA position.
MISSING = {"-", "N", "n", "."}

# --preset -> (iqr_mult, min_bad). Tunes how aggressive the MRate2 fence is.
PRESETS = {
    "strict":    (1.0, 1),
    "default":   (1.5, 2),
    "sensitive": (2.0, 3),
}


# ----------------------------------------------------------------------------
# Small numeric helpers (R quantile type 7, matching qcmsa.R's quantile()).
# ----------------------------------------------------------------------------
def quantile_type7(sorted_vals, p):
    """Linear-interpolation quantile (R type 7 / numpy default) on a sorted list."""
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return float(sorted_vals[0])
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    frac = h - lo
    return sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo])


def upper_fence(values, iqr_mult):
    """Tukey upper fence: Q3 + iqr_mult * (Q3 - Q1)."""
    s = sorted(values)
    q25 = quantile_type7(s, 0.25)
    q75 = quantile_type7(s, 0.75)
    return q75 + iqr_mult * (q75 - q25)


def flag_outliers(counts, min_bad, iqr_mult):
    """Given {key: n_bad}, return (flagged_set, fence) for upper-outlier keys.

    Mirrors qcmsa.R::flag_outliers -- fence is computed on the *non-zero* counts,
    needs >= 4 non-zero values to fire at all, and a key must additionally have
    n_bad >= min_bad so a single bad cell can never trigger removal on its own.
    Returns the fence used (or None when too few non-zero values to fire).
    """
    nonzero = [c for c in counts.values() if c > 0]
    if len(nonzero) < 4:
        return set(), None
    fence = upper_fence(nonzero, iqr_mult)
    return {k for k, c in counts.items() if c > fence and c >= min_bad}, fence


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


_PART_RE = re.compile(r"gene(\d+)\s*=\s*(\d+)\s*-\s*(\d+)")


def parse_partition(path):
    """Return list of (gene_id:int, start:int, end:int), 0-based inclusive."""
    genes = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            m = _PART_RE.search(line)
            if not m:
                continue
            genes.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    return genes


# protal .meta.tsv columns (the header protal writes). We read by header name so
# the extra coverage columns added over time don't shift our indices, and -- the
# bug qcmsa.R had -- the header line is NOT mistaken for data.
META_GENE_COL = "gene_id"
META_SAMPLE_COL = "sample"
META_MRATE2_COL = "multi_rate_vcov2"
META_HCOV_COL = "hcov"               # fraction of gene covered (M3 --gene_min_hcov_frac)
META_DEPTH_COL = "mean_vcov_nonzero"  # mean depth over covered positions (M3 --gene_min_mean_depth)


def load_meta(path, gene_whitelist):
    """Load per-(sample,gene) meta restricted to genes in gene_whitelist.

    Returns (rows, samples_in_order, genes_sorted, cov) where
      rows = list of (sample:str, gene:int, mrate2:float)
      cov  = {(sample, gene): (hcov:float, mean_depth:float)} (empty if the
             coverage columns are absent).
    """
    rows = []
    cov = {}
    with open_maybe_gz(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            si = header.index(META_SAMPLE_COL)
            gi = header.index(META_GENE_COL)
            mi = header.index(META_MRATE2_COL)
        except ValueError as exc:
            raise SystemExit(
                f"qcmsa.py: meta file '{path}' is missing expected column "
                f"({exc}); header was: {header}"
            )
        hi = header.index(META_HCOV_COL) if META_HCOV_COL in header else None
        di = header.index(META_DEPTH_COL) if META_DEPTH_COL in header else None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            gene = int(f[gi])
            if gene not in gene_whitelist:
                continue
            sample = f[si]
            rows.append((sample, gene, float(f[mi])))
            if hi is not None and di is not None:
                try:
                    cov[(sample, gene)] = (float(f[hi]), float(f[di]))
                except ValueError:
                    pass

    samples_seen = []
    seen = set()
    genes_seen = set()
    for sample, gene, _ in rows:
        if sample not in seen:
            seen.add(sample)
            samples_seen.append(sample)
        genes_seen.add(gene)
    return rows, samples_seen, sorted(genes_seen), cov


def coverage_filter(cov, all_genes, hcov_t, depth_t, min_samples):
    """Reproduce protal's M3 coverage gate from the meta coverage columns.

    A (sample,gene) cell passes if hcov >= hcov_t AND mean_depth >= depth_t.
    A gene is dropped if NOT more than min_samples cells pass (protal uses a
    strict '>' on msa_min_samples). Returns (dropped_genes, fail_cells, reason)
    where fail_cells are coverage-failing cells in *surviving* genes (to gap-fill).
    """
    passing = defaultdict(int)
    failing = defaultdict(set)   # gene -> set of failing samples
    for (sample, gene), (h, d) in cov.items():
        if gene not in all_genes:
            continue
        if h >= hcov_t and d >= depth_t:
            passing[gene] += 1
        else:
            failing[gene].add(sample)
    dropped, fail_cells, reason = set(), set(), {}
    for gene in all_genes:
        np = passing.get(gene, 0)
        if np <= min_samples:        # not strictly greater -> dropped (matches protal)
            dropped.add(gene)
            reason[gene] = (np, "M3 coverage: only %d sample(s) pass" % np)
        else:
            for s in failing.get(gene, ()):
                fail_cells.add((s, gene))
    return dropped, fail_cells, reason


def read_fasta(path):
    """Return (names:list, seqs:list) preserving file order."""
    names, seqs = [], []
    cur = []
    with open_maybe_gz(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].strip())
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


def write_fasta(path, names, seqs, width=80):
    with open(path, "w") as fh:
        for name, seq in zip(names, seqs):
            fh.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------------
# Pass 1: MRate2 iterative gene/sample filter (identical to qcmsa.R).
# ----------------------------------------------------------------------------
def mrate2_filter(rows, all_genes, all_samples, min_bad, iqr_mult,
                  gene_abs=0, sample_abs=0, max_iter=100):
    kept_genes = set(all_genes)
    kept_samples = set(all_samples)
    # reason[id] = (n_bad, fence_or_'abs>=N', iteration) recorded when flagged
    gene_reason = {}
    sample_reason = {}

    for it in range(1, max_iter + 1):
        # gene -> number of kept samples with MRate2 > 0
        gene_counts = defaultdict(int)
        for sample, gene, mr in rows:
            if gene in kept_genes and sample in kept_samples and mr > 0:
                gene_counts[gene] += 1
        # genes with zero bad cells still need a (zero) entry for the fence base
        for g in kept_genes:
            gene_counts.setdefault(g, 0)
        tukey_genes, gene_fence = flag_outliers(gene_counts, min_bad, iqr_mult)
        # absolute rule: catch ANY signal above the (often clean-zero) baseline,
        # which the Tukey fence cannot do when the bad items ARE the distribution.
        abs_genes = ({g for g, c in gene_counts.items() if c >= gene_abs}
                     if gene_abs > 0 else set())
        bad_genes = tukey_genes | abs_genes

        # sample -> number of bad genes, computed AFTER removing this round's bad genes
        sample_counts = defaultdict(int)
        for sample, gene, mr in rows:
            if (gene in kept_genes and gene not in bad_genes
                    and sample in kept_samples and mr > 0):
                sample_counts[sample] += 1
        for s in kept_samples:
            sample_counts.setdefault(s, 0)
        tukey_samples, sample_fence = flag_outliers(sample_counts, min_bad, iqr_mult)
        abs_samples = ({s for s, c in sample_counts.items() if c >= sample_abs}
                       if sample_abs > 0 else set())
        bad_samples = tukey_samples | abs_samples

        sys.stderr.write(
            f"  iter {it}: {len(bad_genes)} gene(s) flagged "
            f"({len(tukey_genes)} Tukey, {len(abs_genes - tukey_genes)} abs) | "
            f"{len(bad_samples)} sample(s) flagged "
            f"({len(tukey_samples)} Tukey, {len(abs_samples - tukey_samples)} abs)\n"
        )

        if not bad_genes and not bad_samples:
            break
        for g in bad_genes:
            fence = (gene_fence if g in tukey_genes else f">=abs {gene_abs}")
            gene_reason[g] = (gene_counts[g], fence, it)
        for s in bad_samples:
            fence = (sample_fence if s in tukey_samples else f">=abs {sample_abs}")
            sample_reason[s] = (sample_counts[s], fence, it)
        kept_genes -= bad_genes
        kept_samples -= bad_samples

    filtered_genes = set(all_genes) - kept_genes
    filtered_samples = set(all_samples) - kept_samples
    return (kept_genes, kept_samples, filtered_genes, filtered_samples,
            gene_reason, sample_reason)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def build_argparser():
    p = argparse.ArgumentParser(
        description="QC and filter a strain MSA using protal .meta.tsv metrics."
    )
    p.add_argument("msa", help="MSA FASTA (plain or .gz)")
    p.add_argument("partition", help="RAxML-style partition file")
    p.add_argument("meta", help="protal .meta.tsv (with header)")
    p.add_argument("--prefix", default=None,
                   help="Output prefix (default: MSA path with .fna/.gz stripped)")

    p.add_argument("--preset", choices=sorted(PRESETS),
                   help="Convenience bundle for --iqr-mult/--min-bad "
                        "(strict=1.0/1, default=1.5/2, sensitive=2.0/3). "
                        "Explicit --iqr-mult/--min-bad override the preset.")
    p.add_argument("--iqr-mult", type=float, default=None,
                   help="Tukey IQR multiplier for the MRate2 fence (default 1.5)")
    p.add_argument("--min-bad", type=int, default=None,
                   help="Min bad peers before a gene/sample is removed (default 2)")
    # Absolute multi-allelicity cutoffs -- catch ANY signal above a clean-zero
    # baseline (e.g. conspecific/mixed strains), which the Tukey fence cannot do
    # because the bad items then ARE the distribution. Default 0 = off.
    p.add_argument("--sample-abs-min-bad", type=int, default=0,
                   help="Remove a sample with multi-allelic (MRate2>0) signal in >= this "
                        "many genes, regardless of the Tukey fence. 0=off. Try 1-2 to catch "
                        "mixed/conspecific strains.")
    p.add_argument("--gene-abs-min-bad", type=int, default=0,
                   help="Remove a gene multi-allelic in >= this many samples, regardless of "
                        "the Tukey fence. 0=off.")

    # Coverage gating -- this is where the gene/sample coverage filtering lives
    # (protal emits a raw MSA). Computed from the meta hcov / mean-depth columns.
    # Defaults are ON; set any to 0 to disable that part.
    p.add_argument("--gene-min-hcov", type=float, default=0.3,
                   help="Min fraction of a gene covered for a (sample,gene) cell to "
                        "pass. Default 0.3; 0 disables.")
    p.add_argument("--gene-min-mean-depth", type=float, default=1.0,
                   help="Min mean depth over covered positions for a cell to pass. "
                        "Default 1.0; 0 disables.")
    p.add_argument("--gene-min-samples", type=int, default=3,
                   help="Drop a gene unless MORE than this many samples pass coverage "
                        "(strict >, like protal's old msa_min_samples). Default 3; 0 disables.")
    p.add_argument("--max-mrate2", type=float, default=None,
                   help="Hard per-cell MRate2 cap for the cell-outlier fence "
                        "(default: Tukey fence derived from the data)")

    # M5 step 4b -- site / sequence cleanup
    p.add_argument("--remove-constant", dest="remove_constant",
                   action="store_true", default=True,
                   help="Drop constant sites (default: on)")
    p.add_argument("--keep-constant", dest="remove_constant", action="store_false",
                   help="Keep constant sites")
    p.add_argument("--min-parsimony-samples", type=int, default=2,
                   help="Drop sites where fewer than N samples differ from the "
                        "majority base (default 2; subsumes constant-site removal)")
    p.add_argument("--reapply-hcov", type=int, default=0,
                   help="After gene/site removal, drop sequences with fewer than "
                        "this many valid (non -/N) bases (default 0 = disabled). "
                        "protal passes its --msa_min_hcov here.")
    p.add_argument("--mask-cell-outliers", dest="mask_cells",
                   action="store_true", default=True,
                   help="Mask individual (sample,gene) MRate2 outlier cells with "
                        "'-' (default: on). qcmsa.R only marks these on the plot.")
    p.add_argument("--no-mask-cell-outliers", dest="mask_cells",
                   action="store_false", help="Do not mask cell outliers (qcmsa.R parity)")

    p.add_argument("--plot", action="store_true",
                   help="Also emit a <prefix>.qc.png MRate2 heatmap (needs matplotlib)")
    p.add_argument("--no-summary", dest="summary", action="store_false", default=True,
                   help="Do not write the <prefix>.qcmsa_summary.tsv decision log")
    return p


def resolve_params(args):
    iqr_mult, min_bad = PRESETS["default"]
    if args.preset:
        iqr_mult, min_bad = PRESETS[args.preset]
    if args.iqr_mult is not None:
        iqr_mult = args.iqr_mult
    if args.min_bad is not None:
        min_bad = args.min_bad
    return iqr_mult, min_bad


def main(argv=None):
    args = build_argparser().parse_args(argv)
    iqr_mult, min_bad = resolve_params(args)

    prefix = args.prefix
    if prefix is None:
        # strip ".raw.msa.fna" / ".msa.fna" / ".fna" so the output is <name>.msa.fna
        prefix = re.sub(r"(\.raw)?(\.msa)?\.fna(\.gz)?$", "", args.msa)
    if os.path.abspath(prefix + ".msa.fna") == os.path.abspath(args.msa):
        raise SystemExit("qcmsa.py: output would overwrite the input MSA; pass a "
                         "distinct --prefix (input should be <name>.raw.msa.fna).")

    # --- inputs ---
    partition = parse_partition(args.partition)
    if not partition:
        raise SystemExit(f"qcmsa.py: no genes parsed from partition '{args.partition}'")
    gene_whitelist = {g for g, _, _ in partition}
    sys.stderr.write(f"Partition: {len(partition)} genes\n")

    rows, all_samples, all_genes, cov = load_meta(args.meta, gene_whitelist)
    all_samples = sorted(all_samples)
    sys.stderr.write(
        f"Loaded meta: {len(all_samples)} samples x {len(all_genes)} genes (in MSA)\n"
    )
    sys.stderr.write(f"Params: iqr_mult={iqr_mult} min_bad={min_bad}"
                     + (f" preset={args.preset}" if args.preset else "") + "\n")

    # --- pass 0: optional coverage gate (reproduces protal M3 from meta) ---
    cov_gate = (args.gene_min_hcov > 0 or args.gene_min_mean_depth > 0
                or args.gene_min_samples > 0)
    cov_dropped_genes, cov_fail_cells, cov_reason = set(), set(), {}
    if cov_gate:
        if not cov:
            sys.stderr.write("qcmsa.py: coverage gating requested but meta has no "
                             "hcov/mean_vcov_nonzero columns; skipping coverage gate.\n")
        else:
            cov_dropped_genes, cov_fail_cells, cov_reason = coverage_filter(
                cov, set(all_genes), args.gene_min_hcov,
                args.gene_min_mean_depth, args.gene_min_samples)
            sys.stderr.write(
                f"Coverage gate (hcov>={args.gene_min_hcov}, depth>="
                f"{args.gene_min_mean_depth}, >{args.gene_min_samples} samples): "
                f"dropped {len(cov_dropped_genes)}/{len(all_genes)} genes, "
                f"gap-filled {len(cov_fail_cells)} cell(s)\n")
    # Genes/cells removed by coverage don't participate in the MRate2 stats (they
    # are gaps in the output), so the adaptive fence is computed on covered data.
    cov_survivor_genes = [g for g in all_genes if g not in cov_dropped_genes]
    rows_for_mrate2 = [(s, g, mr) for (s, g, mr) in rows
                       if g not in cov_dropped_genes and (s, g) not in cov_fail_cells]

    # --- pass 1: MRate2 gene/sample filter ---
    (kept_genes, kept_samples, mr_filtered_genes, filtered_samples,
     gene_reason, sample_reason) = mrate2_filter(
        rows_for_mrate2, cov_survivor_genes, all_samples, min_bad, iqr_mult,
        gene_abs=args.gene_abs_min_bad, sample_abs=args.sample_abs_min_bad
    )
    filtered_genes = mr_filtered_genes | cov_dropped_genes
    sys.stderr.write(f"Filtered genes: {len(filtered_genes)} / {len(all_genes)}"
                     f" ({len(cov_dropped_genes)} coverage, {len(mr_filtered_genes)} MRate2)\n")
    sys.stderr.write(f"Filtered samples: {len(filtered_samples)} / {len(all_samples)}\n")

    # --- cell-outlier fence ---
    # NB: qcmsa.R fences over *all* MRate2 values, which collapses to 0 when the
    # data are sparse (mostly zeros) and would then flag every multi-allelic cell.
    # We fence over the non-zero values instead (same rationale as the count
    # fences) and disable masking if the fence is still degenerate.
    nonzero_mrate2 = [mr for _, _, mr in rows if mr > 0]
    if args.max_mrate2 is not None:
        cell_fence = args.max_mrate2
    elif len(nonzero_mrate2) >= 4:
        f = upper_fence(nonzero_mrate2, iqr_mult)
        cell_fence = f if f > 0 else float("inf")
    else:
        cell_fence = float("inf")
    # (sample, gene) cells that exceed the fence and survived the row/col filter
    outlier_cells = {
        (sample, gene)
        for sample, gene, mr in rows
        if mr > cell_fence and sample not in filtered_samples and gene not in filtered_genes
    }
    sys.stderr.write(
        f"Cell fence (MRate2): {cell_fence:.5f} -> {len(outlier_cells)} outlier cell(s)\n"
    )

    # --- read MSA ---
    names, seqs = read_fasta(args.msa)
    if not names:
        raise SystemExit(f"qcmsa.py: empty MSA '{args.msa}'")
    seq_of = dict(zip(names, seqs))
    sys.stderr.write(f"MSA loaded: {len(names)} sequences, {len(seqs[0])} bp\n")

    sample_set = set(all_samples)
    # Keep references (names not in meta) always; drop filtered samples.
    kept_names = [n for n in names if n not in sample_set or n not in filtered_samples]

    # Degenerate MSA (e.g. only the reference survived protal's row filter): nothing
    # meaningful to filter. Warn and skip gracefully so batch/--run_qcmsa runs continue.
    n_sample_seqs = sum(1 for n in kept_names if n in sample_set)
    if n_sample_seqs < 2:
        sys.stderr.write(
            f"qcmsa.py: only {n_sample_seqs} sample sequence(s) in '{args.msa}' "
            "after filtering; nothing to filter - skipping.\n"
        )
        return 0

    # --- surviving genes -> original column ranges (sorted by start) ---
    kept_partition = sorted(
        [(g, s, e) for (g, s, e) in partition if g not in filtered_genes],
        key=lambda t: t[1],
    )
    if not kept_partition:
        sys.stderr.write("qcmsa.py: all genes filtered - no MSA produced (skipping).\n")
        return 0

    # Per surviving column: which gene it belongs to, and its original index.
    col_gene = []
    col_orig = []
    for g, s, e in kept_partition:
        for c in range(s, e + 1):
            col_gene.append(g)
            col_orig.append(c)

    # Subset each kept sequence to the surviving columns; gap-fill coverage-failed
    # cells (M3-equivalent) and, if requested, MRate2 outlier cells.
    msa_rows = []
    for n in kept_names:
        seq = seq_of[n]
        chars = [seq[c] for c in col_orig]
        if n in sample_set:
            for j, g in enumerate(col_gene):
                if (n, g) in cov_fail_cells or (args.mask_cells and (n, g) in outlier_cells):
                    chars[j] = "-"
        msa_rows.append(chars)

    n_cols = len(col_gene)

    # --- reapply-hcov: drop sequences below the valid-base floor AFTER gene removal ---
    # Measured on the gene-filtered alignment (full gene columns), NOT after site
    # cleanup -- otherwise the floor is compared against the tiny informative-only
    # alignment and would drop every sample.
    if args.reapply_hcov > 0:
        keep_idx = []
        for i, (n, row) in enumerate(zip(kept_names, msa_rows)):
            if n in sample_set:  # references are exempt from the floor
                valid = sum(1 for ch in row if ch not in MISSING)
                if valid < args.reapply_hcov:
                    continue
            keep_idx.append(i)
        dropped = len(kept_names) - len(keep_idx)
        if dropped:
            sys.stderr.write(
                f"reapply-hcov={args.reapply_hcov}: dropped {dropped} sequence(s) "
                "below the valid-base floor\n"
            )
        kept_names = [kept_names[i] for i in keep_idx]
        msa_rows = [msa_rows[i] for i in keep_idx]

    # --- pass 2: site cleanup (constant / low-parsimony) ---
    col_keep = [True] * n_cols
    if args.remove_constant or args.min_parsimony_samples > 0:
        for j in range(n_cols):
            counts = Counter()
            for row in msa_rows:
                ch = row[j]
                if ch not in MISSING:
                    counts[ch] += 1
            total = sum(counts.values())
            if total == 0:
                col_keep[j] = False  # all-missing column: nothing to keep
                continue
            majority = max(counts.values())
            minor = total - majority
            distinct = len(counts)
            if args.remove_constant and distinct <= 1:
                col_keep[j] = False
            elif minor < args.min_parsimony_samples:
                col_keep[j] = False

    surviving = [j for j in range(n_cols) if col_keep[j]]
    n_removed_sites = n_cols - len(surviving)
    if not surviving:
        sys.stderr.write(
            "qcmsa.py: all sites removed by cleanup - no MSA produced (skipping). "
            "Consider --keep-constant / --min-parsimony-samples 0.\n"
        )
        return 0

    # Final sequences (column subset).
    final_seqs = ["".join(row[j] for j in surviving) for row in msa_rows]

    sys.stderr.write(
        f"Site cleanup: removed {n_removed_sites} / {n_cols} sites "
        f"({len(surviving)} retained)\n"
    )

    # --- write filtered MSA ---
    msa_out = prefix + ".msa.fna"
    write_fasta(msa_out, kept_names, final_seqs)
    sys.stderr.write(f"Saved: {msa_out}\n")

    # --- write updated partition (recomputed contiguous coordinates) ---
    # Surviving column count per gene, in kept_partition order.
    surv_per_gene = defaultdict(int)
    for j in surviving:
        surv_per_gene[col_gene[j]] += 1
    part_out = prefix + ".partition.txt"
    with open(part_out, "w") as fh:
        new_start = 0
        for g, _, _ in kept_partition:
            length = surv_per_gene.get(g, 0)
            if length == 0:
                continue  # gene lost all its sites in cleanup
            new_end = new_start + length - 1
            fh.write(f"DNA, gene{g} = {new_start}-{new_end}\n")
            new_start = new_end + 1
    sys.stderr.write(f"Saved: {part_out}\n")

    # --- machine-readable decision log (drives the strain_report plots) ---
    if args.summary:
        n_seqs_out = len(kept_names)
        n_ref = sum(1 for n in kept_names if n not in sample_set)
        summary_out = prefix + ".qcmsa_summary.tsv"
        with open(summary_out, "w") as fh:
            fh.write("section\tkey\tvalue\treason\n")
            fh.write(f"param\tiqr_mult\t{iqr_mult}\t\n")
            fh.write(f"param\tmin_bad\t{min_bad}\t\n")
            fh.write(f"param\tpreset\t{args.preset or ''}\t\n")
            fh.write(f"param\tcell_fence\t{cell_fence:.6g}\t\n")
            fh.write(f"count\tsamples_in\t{len(all_samples)}\t\n")
            fh.write(f"count\tsamples_kept\t{len(kept_samples)}\t\n")
            fh.write(f"count\tsamples_filtered\t{len(filtered_samples)}\t\n")
            fh.write(f"count\tgenes_in\t{len(all_genes)}\t\n")
            fh.write(f"count\tgenes_kept\t{len(kept_genes)}\t\n")
            fh.write(f"count\tgenes_filtered\t{len(filtered_genes)}\t\n")
            fh.write(f"count\tgenes_filtered_coverage\t{len(cov_dropped_genes)}\t\n")
            fh.write(f"count\tgenes_filtered_mrate2\t{len(mr_filtered_genes)}\t\n")
            fh.write(f"count\tcoverage_gap_filled_cells\t{len(cov_fail_cells)}\t\n")
            fh.write(f"count\tsites_in\t{n_cols}\t\n")
            fh.write(f"count\tsites_kept\t{len(surviving)}\t\n")
            fh.write(f"count\tsites_removed\t{n_removed_sites}\t\n")
            fh.write(f"count\toutlier_cells\t{len(outlier_cells)}\t\n")
            fh.write(f"count\tseqs_out\t{n_seqs_out}\t\n")
            fh.write(f"count\treference_seqs_out\t{n_ref}\t\n")
            for g in sorted(filtered_genes):
                if g in cov_dropped_genes:
                    nb, txt = cov_reason.get(g, (None, "M3 coverage"))
                    fh.write(f"gene_filtered\t{g}\t{nb}\t{txt}\n")
                    continue
                nb, fence, it = gene_reason.get(g, (None, None, None))
                fs = (f"{fence:.3g}" if isinstance(fence, (int, float)) else fence)
                reason = (f"multi-allelic: {nb} samples > {fs} "
                          f"(iter {it})") if nb is not None else "multi-allelic outlier"
                fh.write(f"gene_filtered\t{g}\t{nb}\t{reason}\n")
            for s in sorted(filtered_samples):
                nb, fence, it = sample_reason.get(s, (None, None, None))
                fs = (f"{fence:.3g}" if isinstance(fence, (int, float)) else fence)
                reason = (f"multi-allelic: {nb} genes > {fs} "
                          f"(iter {it})") if nb is not None else "multi-allelic outlier"
                fh.write(f"sample_filtered\t{s}\t{nb}\t{reason}\n")
            for s, g in sorted(outlier_cells):
                fh.write(f"cell_outlier\t{s}|gene{g}\t\tMRate2 > cell_fence {cell_fence:.3g}\n")
        sys.stderr.write(f"Saved: {summary_out}\n")

    if args.plot:
        try:
            make_plot(prefix, rows, all_genes, all_samples,
                      filtered_genes, filtered_samples, outlier_cells)
        except Exception as exc:  # plotting must never break the filter
            sys.stderr.write(f"qcmsa.py: --plot failed ({exc}); skipping heatmap\n")

    sys.stderr.write("Done.\n")
    return 0


def make_plot(prefix, rows, all_genes, all_samples,
              filtered_genes, filtered_samples, outlier_cells):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gidx = {g: i for i, g in enumerate(all_genes)}
    sidx = {s: i for i, s in enumerate(all_samples)}
    grid = [[float("nan")] * len(all_genes) for _ in all_samples]
    for sample, gene, mr in rows:
        if sample in sidx and gene in gidx:
            grid[sidx[sample]][gidx[gene]] = mr

    fig, ax = plt.subplots(figsize=(max(6, len(all_genes) * 0.06 + 1.5),
                                    max(3, len(all_samples) * 0.06 + 1.5)))
    im = ax.imshow(grid, aspect="auto", cmap="Blues", interpolation="nearest")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample")
    ax.set_title("MRate2 (red label = filtered, x = outlier cell)")
    ax.set_yticks(range(len(all_samples)))
    ax.set_yticklabels(all_samples, fontsize=4)
    for tick, s in zip(ax.get_yticklabels(), all_samples):
        if s in filtered_samples:
            tick.set_color("red")
    ax.set_xticks([])
    for sample, gene in outlier_cells:
        if sample in sidx and gene in gidx:
            ax.plot(gidx[gene], sidx[sample], marker="x", color="red", markersize=2)
    fig.colorbar(im, ax=ax, shrink=0.5, label="MRate2")
    out = prefix + ".qc.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    sys.stderr.write(f"Saved: {out}\n")


if __name__ == "__main__":
    sys.exit(main())
