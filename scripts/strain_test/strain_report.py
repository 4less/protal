#!/usr/bin/env python3
"""strain_report.py - self-contained QC report for a protal strain run + qcmsa.

Reads a protal strain output directory (the MSAs / partitions / meta / snp_stats
that protal writes) and the qcmsa.py outputs produced from it, then writes a
single self-contained **report.html** (plus report.md and summary.tsv):

  * counts samples and genes in each MSA before (protal default) and after qcmsa
  * automated PASS/WARN/FAIL sanity checks ("did something go wrong?")
  * plot (a): which SNPs were filtered out and why (from .snp_stats.tsv, M1)
  * plot (b): gene- and sample-level qcmsa filtering with reasons (.qcmsa_summary.tsv)
  * per-species MRate2 (multi-allelicity) heatmaps with filtered rows/cols flagged

Pure standard library: plots are rendered as inline SVG, so there is no
matplotlib / numpy / pandas dependency and the HTML is fully portable.
"""

import argparse
import glob
import html
import os
import sys
from collections import defaultdict


# ----------------------------------------------------------------------------
# Input readers
# ----------------------------------------------------------------------------
def count_fasta_seqs(path):
    """Return (total_seqs, sample_seqs). Names containing 'reference' are refs."""
    if not os.path.exists(path):
        return None, None
    total, refs = 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                total += 1
                if "reference" in line.lower():
                    refs += 1
    return total, total - refs


def count_partition_genes(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def load_summary(path):
    if not os.path.exists(path):
        return None
    out = {"param": {}, "count": {}, "genes": [], "samples": [], "cells": []}
    for r in read_tsv(path):
        sec = r["section"]
        if sec == "param":
            out["param"][r["key"]] = r["value"]
        elif sec == "count":
            v = r["value"]
            out["count"][r["key"]] = int(v) if v.lstrip("-").isdigit() else v
        elif sec == "gene_filtered":
            out["genes"].append((r["key"], r["value"], r["reason"]))
        elif sec == "sample_filtered":
            out["samples"].append((r["key"], r["value"], r["reason"]))
        elif sec == "cell_outlier":
            out["cells"].append((r["key"], r["reason"]))
    return out


SNP_FIELDS = [
    ("snps_retained", "retained"),
    ("variants_filtered_qual_sum", "filtered: low phred-sum"),
    ("variants_filtered_obs_cov", "filtered: low coverage"),
]
POS_FIELDS = [
    ("positions_below_min_cov", "below min cov"),
    ("positions_no_coverage", "no coverage"),
    ("valid_positions_removed_by_vcov", "removed by vcov"),
]


def aggregate_snp_stats(path):
    if not os.path.exists(path):
        return None, 0
    agg = defaultdict(float)
    n = 0
    for r in read_tsv(path):
        if "reference" in r.get("sample", "").lower():
            continue
        n += 1
        for col, _ in SNP_FIELDS + POS_FIELDS + [("total_variant_positions", "")]:
            try:
                agg[col] += float(r.get(col, 0) or 0)
            except ValueError:
                pass
    return agg, n


def discover_species(strains_dir):
    metas = sorted(glob.glob(os.path.join(strains_dir, "*.meta.tsv")))
    return [os.path.basename(m)[:-len(".meta.tsv")] for m in metas]


def load_db_gene_counts(path):
    """species -> db_markers (count of marker genes in the species' DB genome)."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    for r in read_tsv(path):
        v = r.get("db_markers", "")
        if v.isdigit():
            out[r["species"]] = int(v)
    return out


# Default protal M3 thresholds (see Options.h). Used if not overridden/parsed.
M3_DEFAULTS = dict(hcov=0.50, depth=3.0, min_samples=3)


def parse_m3_params(log_path):
    """Read the M3 thresholds back from protal's ToString() dump in its run log."""
    if not log_path or not os.path.exists(log_path):
        return None
    keys = {"gene min hcov frac:": "hcov", "gene min mean depth:": "depth",
            "msa min samples:": "min_samples"}
    found = {}
    with open(log_path, errors="ignore") as fh:
        for line in fh:
            for k, name in keys.items():
                if k in line:
                    try:
                        found[name] = float(line.split(k, 1)[1].strip().split()[0])
                    except (ValueError, IndexError):
                        pass
    if len(found) < 3:
        return None
    found["min_samples"] = int(found["min_samples"])
    return found


def compute_protal_filtering(meta_path, hcov_t, depth_t, min_samples):
    """Reconstruct protal's M3 gene-coverage funnel from the per-cell meta table.

    A (sample,gene) cell 'passes' if hcov >= hcov_t AND mean_vcov_nonzero >= depth_t.
    A gene enters the MSA if #passing samples > min_samples (protal uses strict >).
    Returns per-gene/per-cell tallies, or None if the meta file is missing.
    """
    if not os.path.exists(meta_path):
        return None
    obs = defaultdict(int)        # gene -> observed cells
    passed = defaultdict(int)     # gene -> passing cells
    cells = dict(total=0, ok=0, fail_h=0, fail_d=0, fail_both=0)
    samples = set()
    for r in read_tsv(meta_path):
        samples.add(r["sample"])
        g = r["gene_id"]
        try:
            h = float(r.get("hcov", 0) or 0)
            d = float(r.get("mean_vcov_nonzero", 0) or 0)
        except ValueError:
            h, d = 0.0, 0.0
        obs[g] += 1
        cells["total"] += 1
        ok_h, ok_d = h >= hcov_t, d >= depth_t
        if ok_h and ok_d:
            passed[g] += 1
            cells["ok"] += 1
        elif not ok_h and not ok_d:
            cells["fail_both"] += 1
        elif not ok_h:
            cells["fail_h"] += 1
        else:
            cells["fail_d"] += 1
    genes_obs = len(obs)
    genes_in = sum(1 for g in obs if passed[g] > min_samples)
    return dict(genes_observed=genes_obs, genes_in_msa=genes_in,
                genes_dropped=genes_obs - genes_in, cells=cells,
                n_samples=len(samples))


# ----------------------------------------------------------------------------
# Tiny SVG chart helpers (no external deps)
# ----------------------------------------------------------------------------
def _esc(s):
    return html.escape(str(s))


def species_panels(items, colors, ylabel):
    """One row per species: a horizontal 100%-stacked 'relative' bar (with the
    species label on the left) and, to its right, a horizontal 'absolute' bar
    scaled to a shared maximum so magnitudes stay comparable across species.

    items: list of (species_name, [(cat, value), ...]) in stacking order.
    colors: {cat: hex} (also defines legend order)."""
    if not items:
        return "<p><em>No data.</em></p>"

    cats = list(colors.keys())
    maxtotal = max((sum(v for _, v in c) for _, c in items), default=1) or 1

    label_w, rel_w, gap, abs_w, rmar = 215, 235, 30, 300, 60
    row_h, bar_h, top, bot = 23, 15, 64, 30
    width = label_w + rel_w + gap + abs_w + rmar
    height = top + len(items) * row_h + bot

    P = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
         f'font-family="sans-serif" font-size="10">']

    # legend (wraps via tspan-free manual x advance)
    lx = 4
    for cat in cats:
        P.append(f'<rect x="{lx}" y="6" width="11" height="11" fill="{colors[cat]}"/>')
        P.append(f'<text x="{lx+15}" y="16" fill="#333">{_esc(cat)}</text>')
        lx += 15 + len(cat) * 6.0 + 16
    if lx > width:  # widen if legend overflows
        width = int(lx) + 10
        P[0] = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
                f'font-family="sans-serif" font-size="10">')

    rel_x0 = label_w
    abs_x0 = label_w + rel_w + gap
    # column headers
    P.append(f'<text x="{rel_x0+rel_w/2:.0f}" y="{top-16}" text-anchor="middle" '
             f'font-weight="bold" fill="#444">relative (%)</text>')
    P.append(f'<text x="{abs_x0+abs_w/2:.0f}" y="{top-16}" text-anchor="middle" '
             f'font-weight="bold" fill="#444">absolute ({_esc(ylabel)})</text>')

    for i, (name, catvals) in enumerate(items):
        y = top + i * row_h
        ty = y + bar_h - 3
        total = sum(v for _, v in catvals)
        # species label (left)
        P.append(f'<text x="{label_w-8}" y="{ty}" text-anchor="end" fill="#222" '
                 f'font-size="10">{_esc(name)}</text>')
        # relative 100%-stacked bar
        x = rel_x0
        if total > 0:
            for cat, v in catvals:
                w = v / total * rel_w
                if w > 0:
                    P.append(f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="{bar_h}" '
                             f'fill="{colors.get(cat, "#888")}">'
                             f'<title>{_esc(name)} {_esc(cat)}: {v:.0f} ({v/total:.1%})</title></rect>')
                    x += w
        else:
            P.append(f'<rect x="{rel_x0}" y="{y}" width="{rel_w}" height="{bar_h}" fill="#eee"/>')
        # absolute bar (shared scale)
        x = abs_x0
        for cat, v in catvals:
            w = v / maxtotal * abs_w
            if w > 0:
                P.append(f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="{bar_h}" '
                         f'fill="{colors.get(cat, "#888")}">'
                         f'<title>{_esc(name)} {_esc(cat)}: {v:.0f}</title></rect>')
                x += w
        P.append(f'<text x="{x+4:.1f}" y="{ty}" fill="#333" font-size="9">{int(round(total))}</text>')

    # axes under the last row
    base_y = top + len(items) * row_h + 4
    for frac, lab in ((0, "0"), (0.5, "50"), (1.0, "100%")):
        xx = rel_x0 + frac * rel_w
        P.append(f'<line x1="{xx:.1f}" y1="{base_y}" x2="{xx:.1f}" y2="{base_y+4}" stroke="#999"/>')
        P.append(f'<text x="{xx:.1f}" y="{base_y+14}" text-anchor="middle" fill="#777" '
                 f'font-size="8">{lab}</text>')
    for frac in (0, 0.5, 1.0):
        xx = abs_x0 + frac * abs_w
        P.append(f'<line x1="{xx:.1f}" y1="{base_y}" x2="{xx:.1f}" y2="{base_y+4}" stroke="#999"/>')
        P.append(f'<text x="{xx:.1f}" y="{base_y+14}" text-anchor="middle" fill="#777" '
                 f'font-size="8">{int(round(maxtotal*frac))}</text>')
    P.append("</svg>")
    return "\n".join(P)


def svg_gradient_legend(rgb, desired, low_lab="low (0)", high_lab="high"):
    """Small inline colour-scale key: white->rgb swatches, with low/high labels and
    a green check on the 'desired' end. Drawn with discrete rects (no SVG ids)."""
    w, n, x0 = 120, 12, 46
    seg = w / n
    P = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{x0+w+150}" height="30" '
         f'font-family="sans-serif" font-size="9">']
    for k in range(n):
        P.append(f'<rect x="{x0+k*seg:.1f}" y="4" width="{seg+0.6:.1f}" height="10" '
                 f'fill="{_heat_color(k/(n-1), 1.0, rgb)}"/>')
    P.append(f'<rect x="{x0}" y="4" width="{w}" height="10" fill="none" stroke="#bbb"/>')
    P.append(f'<text x="{x0-4}" y="13" text-anchor="end" fill="#666">{_esc(low_lab)}</text>')
    P.append(f'<text x="{x0+w+4}" y="13" fill="#666">{_esc(high_lab)}</text>')
    dx = (x0 if desired == "low" else x0 + w)
    P.append(f'<text x="{dx}" y="27" text-anchor="middle" fill="#2ca02c" '
             f'font-weight="bold">&#10003; desired</text>')
    P.append('<rect x="0" y="2" width="40" height="13" fill="#eeeeee" stroke="#bbb"/>')
    P.append('<text x="44" y="27" fill="#999">grey = gene absent in sample</text>')
    P.append("</svg>")
    return "".join(P)


def _heat_color(v, vmax, rgb=(214, 39, 40)):
    """White -> rgb interpolation; None renders as light grey (absent cell)."""
    if v is None:
        return "#eeeeee"
    if vmax <= 0:
        return "#ffffff"
    t = min(1.0, max(0.0, v / vmax))
    r = int(255 + (rgb[0] - 255) * t)
    g = int(255 + (rgb[1] - 255) * t)
    b = int(255 + (rgb[2] - 255) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def svg_heatmap(samples, genes, cell, filt_samples, filt_genes, title, vmax,
                rgb=(214, 39, 40), fail_below=None, fail_var=None, fmt="%.3f"):
    """Sample x gene heatmap. Optionally outline cells whose value is below a
    filtering threshold (fail_below) to show what drives gap-filling."""
    rows, cols = len(samples), len(genes)
    if rows == 0 or cols == 0:
        return f"<p><em>No heatmap for {_esc(title)}</em></p>"
    cw = max(4, min(14, 900 // cols))
    ch = max(4, min(14, 600 // rows))
    pad_l, pad_t, pad_r, pad_b = 92, 28, 30, 18
    width = pad_l + cols * cw + pad_r
    height = pad_t + rows * ch + pad_b
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'font-family="sans-serif" font-size="9">']
    parts.append(f'<text x="{pad_l}" y="16" font-size="12" font-weight="bold">{_esc(title)}</text>')
    for i, s in enumerate(samples):
        y = pad_t + i * ch
        col = "red" if s in filt_samples else "#333"
        if ch >= 7:
            parts.append(f'<text x="{pad_l-4}" y="{y+ch-1}" text-anchor="end" fill="{col}">{_esc(s)}</text>')
        for j, g in enumerate(genes):
            v = cell.get((s, g))
            fail = (fail_below is not None and v is not None and v < fail_below)
            stroke = ('stroke="#000" stroke-width="0.6"' if fail else 'stroke="#fff" stroke-width="0.2"')
            vtxt = (fmt % v) if v is not None else "NA"
            parts.append(f'<rect x="{pad_l+j*cw}" y="{y}" width="{cw}" height="{ch}" '
                         f'fill="{_heat_color(v, vmax, rgb)}" {stroke}>'
                         f'<title>{_esc(s)} / gene{_esc(g)}: {vtxt}'
                         f'{" (below threshold)" if fail else ""}</title></rect>')
    for j, g in enumerate(genes):
        if g in filt_genes:
            x = pad_l + j * cw + cw / 2
            parts.append(f'<polygon points="{x-3},{pad_t+rows*ch+2} {x+3},{pad_t+rows*ch+2} {x},{pad_t+rows*ch+8}" fill="red"/>')
    parts.append("</svg>")
    return "\n".join(parts)


# Per-(sample,gene) meta variables to render as heatmaps. Fields:
#   key, label, rgb, m3_drive, caption, desired_end, good_meaning, bad_meaning
HEATMAP_VARS = [
    ("hcov", "hcov (fraction of gene covered)", (31, 119, 180), "m3_hcov",
     "drives the M3 gene filter; black-outlined cells fall below --gene_min_hcov_frac and are gap-filled",
     "high", "dark = well-covered gene (kept)", "pale / black-outlined = sparse coverage, dropped by M3"),
    ("mean_vcov_nonzero", "mean depth over covered positions", (44, 160, 44), "m3_depth",
     "drives the M3 gene filter; black-outlined cells fall below --gene_min_mean_depth",
     "high", "dark = deep, confident coverage (kept)", "pale / black-outlined = shallow, dropped by M3"),
    ("vertical_coverage", "vertical coverage (mean depth over whole gene)", (148, 103, 189), None,
     "overall per-gene sequencing depth (context)",
     "high", "dark = more sequencing depth", "pale = barely sequenced"),
    ("multi_rate_vcov2", "MRate2 (multi-allelicity rate)", (214, 39, 40), None,
     "drives qcmsa gene/sample filtering",
     "low", "pale = clean single strain", "dark = mixed strains / contamination, removed by qcmsa"),
    ("filtered_rate_vcov2", "filtered-SNP rate (positions removed by M1 SNP gates)", (255, 127, 14), None,
     "fraction of variant positions M1 discarded",
     "low", "pale = few variants needed filtering", "dark = many low-quality variants discarded"),
]


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strains", required=True, help="protal strain output dir")
    ap.add_argument("--qcmsa", default=None, help="dir with qcmsa outputs (default: --strains)")
    ap.add_argument("--out", required=True, help="report output dir")
    ap.add_argument("--protal-log", default=None,
                    help="protal run log to read M3 thresholds from "
                         "(default: <strains>/../protal_run.log)")
    ap.add_argument("--gene-min-hcov-frac", type=float, default=None)
    ap.add_argument("--gene-min-mean-depth", type=float, default=None)
    ap.add_argument("--msa-min-samples", type=int, default=None)
    ap.add_argument("--db-gene-counts", default=None,
                    help="db_gene_counts.tsv (default: <strains>/../db_gene_counts.tsv)")
    args = ap.parse_args(argv)
    qcmsa_dir = args.qcmsa or args.strains
    os.makedirs(args.out, exist_ok=True)

    db_counts_path = args.db_gene_counts or os.path.join(args.strains, os.pardir, "db_gene_counts.tsv")
    db_markers = load_db_gene_counts(db_counts_path)

    # Resolve M3 thresholds: explicit CLI > parsed from protal log > defaults.
    log_path = args.protal_log or os.path.join(args.strains, os.pardir, "protal_run.log")
    m3 = parse_m3_params(log_path) or dict(M3_DEFAULTS)
    if args.gene_min_hcov_frac is not None:
        m3["hcov"] = args.gene_min_hcov_frac
    if args.gene_min_mean_depth is not None:
        m3["depth"] = args.gene_min_mean_depth
    if args.msa_min_samples is not None:
        m3["min_samples"] = args.msa_min_samples

    species = discover_species(args.strains)
    rows, checks = [], []
    pf_by_species = {}  # species -> protal M3 filtering tallies

    def check(level, sp, msg):
        checks.append((level, sp, msg))

    heat_data = []  # (sp, cell, samples, genes, filt_samples, filt_genes, vmax)

    for sp in species:
        meta = os.path.join(args.strains, sp + ".meta.tsv")
        base_part = os.path.join(args.strains, sp + ".partition.txt")
        pergene = os.path.join(args.strains, sp + ".pergene_filtered.msa.fna")
        base_msa = os.path.join(args.strains, sp + ".msa.fna")
        before_msa = pergene if os.path.exists(pergene) else base_msa

        meta_samples, meta_genes = set(), set()
        var_cells = {vk: {} for vk, *_ in HEATMAP_VARS}   # var -> {(sample,gene): value}
        var_vmax = {vk: 0.0 for vk, *_ in HEATMAP_VARS}
        for r in read_tsv(meta):
            s, g = r["sample"], r["gene_id"]
            meta_samples.add(s); meta_genes.add(g)
            for vk, *_ in HEATMAP_VARS:
                try:
                    v = float(r.get(vk, 0) or 0)
                except ValueError:
                    v = 0.0
                var_cells[vk][(s, g)] = v
                var_vmax[vk] = max(var_vmax[vk], v)
        cell = var_cells["multi_rate_vcov2"]      # MRate2 map (back-compat)
        vmax = var_vmax["multi_rate_vcov2"]

        b_seqs, b_samp = count_fasta_seqs(before_msa)
        b_genes = count_partition_genes(base_part)
        summ = load_summary(os.path.join(qcmsa_dir, sp + ".qcmsa_summary.tsv"))
        a_seqs, a_samp = count_fasta_seqs(os.path.join(qcmsa_dir, sp + ".filtered.msa.fna"))
        a_genes = count_partition_genes(os.path.join(qcmsa_dir, sp + ".filtered.partition.txt"))
        snp_agg, snp_n = aggregate_snp_stats(os.path.join(args.strains, sp + ".snp_stats.tsv"))

        pf = compute_protal_filtering(meta, m3["hcov"], m3["depth"], m3["min_samples"])
        pf_by_species[sp] = pf

        row = dict(species=sp, meta_samples=len(meta_samples), meta_genes=len(meta_genes),
                   before_seqs=b_seqs, before_sample_seqs=b_samp, before_genes=b_genes,
                   after_seqs=a_seqs, after_sample_seqs=a_samp, after_genes=a_genes,
                   qcmsa_ran=summ is not None)
        row["db_markers"] = db_markers.get(sp)
        if pf:
            row["genes_observed"] = pf["genes_observed"]
            row["genes_in_msa"] = pf["genes_in_msa"]
            row["protal_genes_dropped"] = pf["genes_dropped"]
            # sample sequences dropped by protal's per-sequence hcov floor (msa_min_hcov)
            row["protal_samples_dropped"] = (pf["n_samples"] - b_samp
                                             if b_samp is not None else None)
        if summ:
            for k in ("sites_in", "sites_kept", "genes_filtered", "samples_filtered", "outlier_cells"):
                row[k] = summ["count"].get(k)
        rows.append(row)

        if b_seqs and len(meta_samples) >= 2 and len(meta_genes) >= 2:
            fs = {s for s, _, _ in summ["samples"]} if summ else set()
            fg = {g for g, _, _ in summ["genes"]} if summ else set()
            heat_data.append(dict(
                species=sp, samples=sorted(meta_samples),
                genes=sorted(meta_genes, key=lambda g: int(g) if g.isdigit() else g),
                filt_samples=fs, filt_genes=fg,
                var_cells=var_cells, var_vmax=var_vmax))

        # ---- checks ----
        if b_seqs is None:
            check("WARN", sp, "no MSA produced by protal (too few samples/genes for strain reconstruction)")
            continue
        if b_samp is not None and b_samp < 2:
            check("WARN", sp, f"protal MSA has only {b_samp} sample sequence(s); not usable for trees")
        elif b_seqs < 4:
            check("WARN", sp, f"protal MSA has {b_seqs} sequences (<4); too few for a meaningful tree")
        if snp_agg and snp_agg.get("total_variant_positions", 0) > 0 and snp_agg.get("snps_retained", 0) == 0:
            check("WARN", sp, f"ALL {int(snp_agg['total_variant_positions'])} variant positions filtered by "
                              f"SNP filters (M1) across {snp_n} samples - no SNPs retained")
        dbm = db_markers.get(sp)
        if pf and dbm and dbm > 0:
            unhit = dbm - pf["genes_observed"]
            if unhit / dbm >= 0.4:
                check("WARN", sp,
                      f"only {pf['genes_observed']}/{dbm} marker genes in the DB genome got any "
                      f"read hits ({unhit} unhit, {unhit/dbm:.0%}) - abundance/coverage-limited, "
                      "before M3 even applies")
        if pf and pf["genes_observed"] > 0:
            frac = pf["genes_dropped"] / pf["genes_observed"]
            if frac >= 0.4:
                check("WARN", sp,
                      f"protal M3 gene-coverage filter dropped {pf['genes_dropped']}/{pf['genes_observed']} "
                      f"observed genes ({frac:.0%}); low-coverage species - consider lowering "
                      f"--gene_min_hcov_frac ({m3['hcov']}) / --gene_min_mean_depth ({m3['depth']}) "
                      f"or --msa_min_samples ({m3['min_samples']})")
        if summ is None:
            if b_samp and b_samp >= 2:
                check("WARN", sp, "qcmsa produced no output despite >=2 sample sequences")
        else:
            if a_genes is not None and b_genes is not None and a_genes > b_genes:
                check("FAIL", sp, f"qcmsa kept MORE genes ({a_genes}) than input ({b_genes})")
            if a_samp is not None and b_samp is not None and a_samp > b_samp:
                check("FAIL", sp, f"qcmsa kept MORE sample seqs ({a_samp}) than input ({b_samp})")
            si, sk = summ["count"].get("sites_in"), summ["count"].get("sites_kept")
            if isinstance(si, int) and isinstance(sk, int) and si > 0 and sk / si < 0.01:
                check("WARN", sp, f"qcmsa site cleanup kept only {sk}/{si} sites ({sk/si:.1%}); "
                                  "very aggressive (mostly invariant marker positions)")
            if isinstance(a_seqs, int) and a_seqs < 4:
                check("WARN", sp, f"after qcmsa only {a_seqs} sequences remain (<4)")

    write_tables(args.out, rows)
    write_markdown(args.out, rows, checks, species)
    write_html(args.out, rows, checks, species, args.strains, qcmsa_dir, heat_data,
               pf_by_species, m3, db_markers)

    n_fail = sum(1 for l, _, _ in checks if l == "FAIL")
    n_warn = sum(1 for l, _, _ in checks if l == "WARN")
    sys.stderr.write(f"Report written to {os.path.join(args.out, 'report.html')} "
                     f"({n_fail} FAIL, {n_warn} WARN)\n")
    return 1 if n_fail else 0


def write_tables(out, rows):
    cols = ["species", "meta_samples", "meta_genes", "before_seqs", "before_sample_seqs",
            "before_genes", "after_seqs", "after_sample_seqs", "after_genes",
            "sites_in", "sites_kept", "genes_filtered", "samples_filtered", "outlier_cells", "qcmsa_ran"]
    with open(os.path.join(out, "summary.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def write_markdown(out, rows, checks, species):
    n_msa = sum(1 for r in rows if r["before_seqs"])
    L = ["# Protal strain run + qcmsa QC report", "",
         f"- Species with meta table: **{len(species)}**",
         f"- Species with a protal MSA: **{n_msa}**",
         f"- Species with qcmsa output: **{sum(1 for r in rows if r['qcmsa_ran'])}**", ""]
    nf = sum(1 for l, _, _ in checks if l == "FAIL")
    nw = sum(1 for l, _, _ in checks if l == "WARN")
    L.append(f"## Automated checks: {nf} FAIL, {nw} WARN")
    for lvl in ("FAIL", "WARN"):
        for l, sp, msg in checks:
            if l == lvl:
                L.append(f"- **{lvl}** [{sp}]: {msg}")
    if not checks:
        L.append("All checks passed.")
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def write_html(out, rows, checks, species, strains_dir, qcmsa_dir, heat_data,
               pf_by_species=None, m3=None, db_markers=None):
    pf_by_species = pf_by_species or {}
    m3 = m3 or dict(M3_DEFAULTS)
    db_markers = db_markers or {}

    # ---- protal-internal gene/sample filtering (M3) -- per-species panels ----
    pf_items = [(sp, pf) for sp in species if (pf := pf_by_species.get(sp))]
    pf_items.sort(key=lambda t: -(t[1]["genes_observed"]))
    gene_funnel_colors = {"genes in MSA": "#2ca02c", "genes dropped (M3 coverage)": "#d62728",
                          "no read hits (not observed)": "#bbbbbb"}

    def funnel_cats(sp, pf):
        cats = [("genes in MSA", pf["genes_in_msa"]),
                ("genes dropped (M3 coverage)", pf["genes_dropped"])]
        dbm = db_markers.get(sp)
        if dbm:
            cats.append(("no read hits (not observed)", max(0, dbm - pf["genes_observed"])))
        return cats

    funnel_ylabel = ("marker genes in DB genome" if db_markers else "observed genes")
    svg_c_genes = species_panels(
        [(sp.replace("s__", ""), funnel_cats(sp, pf)) for sp, pf in pf_items],
        gene_funnel_colors, funnel_ylabel)
    cell_colors = {"pass coverage": "#2ca02c", "fail: low hcov": "#ff7f0e",
                   "fail: low depth": "#1f77b4", "fail: low hcov+depth": "#d62728"}
    svg_c_cells = species_panels(
        [(sp.replace("s__", ""),
          [("pass coverage", pf["cells"]["ok"]),
           ("fail: low hcov", pf["cells"]["fail_h"]),
           ("fail: low depth", pf["cells"]["fail_d"]),
           ("fail: low hcov+depth", pf["cells"]["fail_both"])])
         for sp, pf in pf_items], cell_colors, "(sample,gene) cells")

    # ---- plot (a): SNP filtering -- one auto-scaled panel per species ----
    snp_rows = []
    for sp in species:
        agg, _ = aggregate_snp_stats(os.path.join(strains_dir, sp + ".snp_stats.tsv"))
        if agg and agg.get("total_variant_positions", 0) > 0:
            snp_rows.append((sp, agg))
    snp_rows.sort(key=lambda d: -d[1].get("total_variant_positions", 0))
    a_colors = {"retained": "#2ca02c", "filtered: low phred-sum": "#ff7f0e",
                "filtered: low coverage": "#d62728"}
    a_items = [(sp.replace("s__", ""), [(lab, agg.get(col, 0)) for col, lab in SNP_FIELDS])
               for sp, agg in snp_rows]
    svg_a = species_panels(a_items, a_colors, "variant positions (Σ samples)")
    p_colors = {"below min cov": "#1f77b4", "no coverage": "#9467bd", "removed by vcov": "#8c564b"}
    p_items = [(sp.replace("s__", ""), [(lab, agg.get(col, 0)) for col, lab in POS_FIELDS])
               for sp, agg in snp_rows]
    svg_a2 = species_panels(p_items, p_colors, "positions (Σ samples)")

    # ---- plot (b): gene & sample filtering -- one auto-scaled panel per species ----
    b_rows = []
    for sp in species:
        s = load_summary(os.path.join(qcmsa_dir, sp + ".qcmsa_summary.tsv"))
        if s:
            b_rows.append((sp, s["count"]))
    g_colors = {"genes kept": "#2ca02c", "genes filtered (MRate2 outlier)": "#d62728"}
    svg_b_genes = species_panels(
        [(sp.replace("s__", ""), [("genes kept", c.get("genes_kept", 0)),
                                  ("genes filtered (MRate2 outlier)", c.get("genes_filtered", 0))])
         for sp, c in b_rows], g_colors, "genes")
    s_colors = {"samples kept": "#1f77b4", "samples filtered (MRate2 outlier)": "#d62728"}
    svg_b_samples = species_panels(
        [(sp.replace("s__", ""), [("samples kept", c.get("samples_kept", 0)),
                                  ("samples filtered (MRate2 outlier)", c.get("samples_filtered", 0))])
         for sp, c in b_rows], s_colors, "samples")

    # ---- per-species filtering-variable heatmaps (one block per species) ----
    heat_blocks = []
    for hd in heat_data:
        s, g = hd["samples"], hd["genes"]
        fs, fg = hd["filt_samples"], hd["filt_genes"]
        svgs = []
        for vk, label, rgb, drive, caption, desired, good, bad in HEATMAP_VARS:
            vmaxv = hd["var_vmax"][vk]
            fail_below = m3["hcov"] if drive == "m3_hcov" else (
                m3["depth"] if drive == "m3_depth" else None)
            fmt = "%.3f" if vmaxv <= 5 else "%.0f"
            high_lab = f"high ({vmaxv:.3g})"
            svgs.append(
                f'<div style="margin:8px 0"><div style="font-size:12px">'
                f'<b>{_esc(label)}</b> <span class="muted">&mdash; {_esc(caption)}</span></div>'
                f'<div class="muted" style="font-size:11px">'
                f'<span style="color:#2ca02c">&#10003; {_esc(good)}</span> &nbsp;&middot;&nbsp; '
                f'<span style="color:#d62728">&#10007; {_esc(bad)}</span></div>'
                + svg_gradient_legend(rgb, desired, high_lab=high_lab)
                + f'<div class="grid">'
                + svg_heatmap(s, g, hd["var_cells"][vk], fs, fg, "",
                              max(1e-9, vmaxv), rgb=rgb, fail_below=fail_below, fmt=fmt)
                + "</div></div>")
        heat_blocks.append(f'<h3>{_esc(hd["species"].replace("s__", ""))} '
                           f'({len(s)} samples &times; {len(g)} genes)</h3>' + "".join(svgs))

    # ---- collect reasons for the filtering tables ----
    reason_blocks = []
    for sp in species:
        s = load_summary(os.path.join(qcmsa_dir, sp + ".qcmsa_summary.tsv"))
        if not s or (not s["genes"] and not s["samples"]):
            continue
        items = []
        for g, nb, reason in s["genes"]:
            items.append(f"<tr><td>gene {_esc(g)}</td><td>{_esc(reason)}</td></tr>")
        for sm, nb, reason in s["samples"]:
            items.append(f"<tr><td>sample {_esc(sm)}</td><td>{_esc(reason)}</td></tr>")
        reason_blocks.append(f"<h4>{_esc(sp)}</h4><table><tr><th>filtered</th><th>reason</th></tr>"
                             + "".join(items) + "</table>")

    nf = sum(1 for l, _, _ in checks if l == "FAIL")
    nw = sum(1 for l, _, _ in checks if l == "WARN")

    def check_rows():
        if not checks:
            return '<p class="ok">All checks passed.</p>'
        out = []
        for lvl in ("FAIL", "WARN"):
            for l, sp, msg in checks:
                if l == lvl:
                    out.append(f'<div class="chk {l.lower()}"><b>{l}</b> '
                               f'[{_esc(sp)}]: {_esc(msg)}</div>')
        return "\n".join(out)

    # summary table
    trs = []
    for r in sorted(rows, key=lambda x: -(x["before_seqs"] or 0)):
        if not r["before_seqs"]:
            continue
        sites = (f'{r.get("sites_in","?")}&rarr;{r.get("sites_kept","?")}' if r["qcmsa_ran"] else "—")
        obs = r.get("genes_observed")
        dbm = r.get("db_markers")
        hit_pct = (f"{obs/dbm:.0%}" if (dbm and obs is not None) else "—")
        kept_pct = (f"{r['before_genes']/obs:.0%}" if obs else "—")
        trs.append(
            f"<tr><td>{_esc(r['species'])}</td><td>{r['meta_samples']}</td>"
            f"<td>{dbm if dbm is not None else '—'}</td>"
            f"<td>{r.get('genes_observed','—')}</td><td>{hit_pct}</td>"
            f"<td>{r['before_genes']}</td><td>{kept_pct}</td><td>{r['before_seqs']}</td>"
            f"<td>{r.get('after_genes','—')}</td><td>{r.get('after_seqs','—')}</td>"
            f"<td>{sites}</td><td>{r.get('genes_filtered','—')}</td>"
            f"<td>{r.get('samples_filtered','—')}</td></tr>")
    no_msa = [r["species"] for r in rows if not r["before_seqs"]]

    n_msa = sum(1 for r in rows if r["before_seqs"])
    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Protal strain + qcmsa QC report</title>
<style>
 body{{font-family:sans-serif;max-width:1100px;margin:24px auto;color:#222;padding:0 16px}}
 h1{{border-bottom:2px solid #444}} h2{{margin-top:32px;border-bottom:1px solid #ccc}}
 table{{border-collapse:collapse;margin:8px 0;font-size:13px}}
 td,th{{border:1px solid #ccc;padding:3px 8px;text-align:right}} th{{background:#f0f0f0}}
 td:first-child,th:first-child{{text-align:left}}
 .chk{{padding:5px 9px;margin:3px 0;border-radius:4px}}
 .fail{{background:#fdd;border-left:4px solid #d62728}}
 .warn{{background:#fff6db;border-left:4px solid #ff7f0e}}
 .ok{{color:#2ca02c}} .muted{{color:#777}} svg{{max-width:100%;height:auto}}
 .grid{{overflow-x:auto}}
</style></head><body>
<h1>Protal strain run + qcmsa QC report</h1>
<p>Species with meta table: <b>{len(species)}</b> &nbsp;|&nbsp;
   with a protal MSA: <b>{n_msa}</b> &nbsp;|&nbsp;
   with qcmsa output: <b>{sum(1 for r in rows if r['qcmsa_ran'])}</b></p>

<h2>Automated checks &mdash; {nf} FAIL, {nw} WARN</h2>
{check_rows()}

<h2>Genes &amp; samples through the pipeline</h2>
<p class="muted">Two filtering stages. <b>protal</b> (M1&ndash;M4) decides which genes/samples enter
the MSA at all (columns up to "protal seqs"); <b>qcmsa</b> (M5) then post-filters that MSA
(columns from "qcmsa genes"). "genes observed" = genes seen in &ge;1 sample (the candidate set);
"protal genes" = genes that passed the M3 coverage filter and entered the MSA.</p>
<table>
<tr><th rowspan="2">species</th><th rowspan="2">meta samples</th>
<th colspan="5">protal (M1&ndash;M4)</th><th rowspan="2">protal seqs</th>
<th colspan="2">qcmsa (M5)</th><th rowspan="2">sites in&rarr;kept</th>
<th rowspan="2">qcmsa genes filt</th><th rowspan="2">qcmsa samples filt</th></tr>
<tr><th>DB markers</th><th>genes observed</th><th>% markers hit</th>
<th>genes in MSA</th><th>% kept (M3)</th>
<th>genes</th><th>seqs</th></tr>
{''.join(trs)}
</table>
<p class="muted">"DB markers" = marker genes present in the species' reference genome (the true
denominator). "genes observed" = got &ge;1 read hit (the rest are abundance-limited, lost
before M3). "% kept (M3)" = of observed genes, the fraction that passed the M3 coverage gate.</p>
<p class="muted">No protal MSA (too few samples/genes to reconstruct a strain): {_esc(', '.join(no_msa) or '(none)')}</p>

<h2>(a) Which SNPs were filtered out, and why &mdash; protal M1</h2>
<p class="muted">Per species, summed across samples. Variants are filtered by protal's M1 SNP gates:
low cumulative phred-sum (<code>--snp_min_phred_sum</code>) or insufficient supporting reads
(<code>--snp_min_cov</code>). The second panel shows why positions had no callable variant.</p>
<p class="muted"><b>Colours:</b>
<span style="color:#2ca02c">&#10003; green = retained</span> (desired &mdash; variants kept for the tree);
<span style="color:#ff7f0e">orange</span> / <span style="color:#d62728">red = filtered out</span>
(undesired loss, but correct when the evidence is weak). A tall green bar is good.</p>
<h3>Variant fate (retained vs filtered by M1 gates)</h3>
<div class="grid">{svg_a}</div>
<h3>Why positions had no callable variant</h3>
<div class="grid">{svg_a2}</div>

<h2>(b) Gene &amp; sample filtering inside protal &mdash; M3 coverage</h2>
<p class="muted">A (sample,gene) cell enters the MSA only if it passes the M3 coverage gate
(<code>--gene_min_hcov_frac</code> = {m3['hcov']}: fraction of gene covered &ge;1 read, AND
<code>--gene_min_mean_depth</code> = {m3['depth']}: mean depth over covered positions). A gene
enters the MSA only if more than <code>--msa_min_samples</code> = {m3['min_samples']} samples pass.
Per species, auto-scaled to its own totals.</p>
<p class="muted"><b>Colours:</b>
<span style="color:#2ca02c">&#10003; green = in MSA</span> (desired);
<span style="color:#d62728">&#10007; red = dropped by M3</span> (observed but coverage too low);
<span style="color:#999">grey = no read hits</span> (never observed &mdash; abundance-limited, the worst case).
More green = better. For the coverage cells, green = passes; orange/blue/red = fails (the colour says
<i>why</i>) and those cells get gap-filled.</p>
<h3>Gene funnel: DB markers &rarr; observed &rarr; kept vs dropped by M3</h3>
<div class="grid">{svg_c_genes}</div>
<h3>Per-(sample,gene) coverage cells: pass vs why they failed</h3>
<div class="grid">{svg_c_cells}</div>

<h2>(c) Gene &amp; sample filtering (qcmsa, milestone M5)</h2>
<p class="muted">qcmsa removes genes/samples that are multi-allelicity (MRate2) outliers via the
iterative Tukey-IQR rule. <b>Colours:</b>
<span style="color:#2ca02c">green</span>/<span style="color:#1f77b4">blue = kept</span> (desired);
<span style="color:#d62728">&#10007; red = removed as a contamination/mixed-strain outlier</span>.
On clean data most bars are fully green/blue (little to remove).</p>
<h3>Gene filtering</h3>
<div class="grid">{svg_b_genes}</div>
<h3>Sample filtering</h3>
<div class="grid">{svg_b_samples}</div>

<h2>(d) Filtering variables per species (sample &times; gene heatmaps)</h2>
<p class="muted">For each species, the per-(sample,gene) variables that drive filtering.
Rows = samples, columns = observed genes (light grey = gene absent in that sample).
Red row labels / red &#9650; columns mark samples / genes removed by qcmsa. Black cell
outlines mark cells below the M3 coverage thresholds
(hcov &lt; {m3['hcov']} or mean depth &lt; {m3['depth']}) that get gap-filled.
<b>Each panel has its own colour key</b> with a &#10003; on the desired end &mdash; note the
direction flips: for the coverage variables <i>darker = better</i>, but for MRate2 and the
filtered-SNP rate <i>darker = worse</i>.</p>
{''.join(heat_blocks)}

<h2>Filtering reasons (detail)</h2>
{''.join(reason_blocks) or '<p class="muted">No genes or samples were filtered by qcmsa on this run.</p>'}

<p class="muted">Generated by scripts/strain_test/strain_report.py</p>
</body></html>"""
    with open(os.path.join(out, "report.html"), "w") as fh:
        fh.write(doc)


if __name__ == "__main__":
    sys.exit(main())
