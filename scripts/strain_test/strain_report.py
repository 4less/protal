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


# ----------------------------------------------------------------------------
# Tiny SVG chart helpers (no external deps)
# ----------------------------------------------------------------------------
def _esc(s):
    return html.escape(str(s))


def svg_species_stack(name, cats, colors, ylabel, width=250, plot_h=230):
    """One vertical stacked bar for a single species, auto-scaled to its own total.

    cats: list of (category_name, value) in stacking order (bottom first)."""
    total = sum(v for _, v in cats)
    vmax = total or 1
    pad_t, pad_b, pad_l, pad_r = 26, 64, 56, 12
    height = pad_t + plot_h + pad_b
    bw = 78
    bx = pad_l + 24
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'font-family="sans-serif" font-size="11">']
    parts.append(f'<text x="{width/2:.0f}" y="16" text-anchor="middle" font-size="11" '
                 f'font-weight="bold">{_esc(name)}</text>')
    # y axis with its OWN scale (0 .. species total)
    parts.append(f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t+plot_h}" stroke="#888"/>')
    for frac in (0, 0.5, 1.0):
        y = pad_t + plot_h * (1 - frac)
        val = int(round(vmax * frac))
        parts.append(f'<line x1="{pad_l-4}" y1="{y:.1f}" x2="{pad_l}" y2="{y:.1f}" stroke="#888"/>')
        parts.append(f'<text x="{pad_l-7}" y="{y+4:.1f}" text-anchor="end" fill="#555" '
                     f'font-size="9">{val}</text>')
    parts.append(f'<text x="13" y="{pad_t+plot_h/2:.0f}" text-anchor="middle" fill="#555" '
                 f'font-size="9" transform="rotate(-90 13 {pad_t+plot_h/2:.0f})">{_esc(ylabel)}</text>')
    y0 = pad_t + plot_h
    for lab, v in cats:
        h = (v / vmax) * plot_h if vmax else 0
        if h > 0:
            parts.append(f'<rect x="{bx}" y="{y0-h:.1f}" width="{bw}" height="{h:.1f}" '
                         f'fill="{colors.get(lab, "#888")}">'
                         f'<title>{_esc(lab)}: {v:.0f} ({v/vmax:.1%})</title></rect>')
        y0 -= h
    parts.append(f'<text x="{bx+bw/2:.0f}" y="{pad_t-4}" text-anchor="middle" fill="#333" '
                 f'font-size="9">Σ={int(round(total))}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def species_panels(items, colors, ylabel):
    """items: list of (species_name, [(cat, value), ...]). Renders one auto-scaled
    mini stacked bar per species in a flex row, with a shared legend on top."""
    if not items:
        return "<p><em>No data.</em></p>"
    legend = "".join(
        f'<span style="white-space:nowrap"><span style="display:inline-block;width:11px;'
        f'height:11px;background:{c};vertical-align:middle;margin:0 4px 0 12px"></span>'
        f'{_esc(n)}</span>' for n, c in colors.items())
    svgs = "".join('<div style="margin:2px 6px">' + svg_species_stack(n, cats, colors, ylabel)
                   + "</div>" for n, cats in items)
    return (f'<div style="margin:4px 0">{legend}</div>'
            f'<div style="display:flex;flex-wrap:wrap;align-items:flex-end">{svgs}</div>')


def _heat_color(v, vmax):
    if v is None:
        return "#eeeeee"
    if vmax <= 0:
        return "#ffffff"
    t = min(1.0, v / vmax)
    # white -> red
    r = 255
    g = int(255 * (1 - t))
    b = int(255 * (1 - t))
    return f"#{r:02x}{g:02x}{b:02x}"


def svg_heatmap(samples, genes, cell, filt_samples, filt_genes, title, vmax):
    rows, cols = len(samples), len(genes)
    if rows == 0 or cols == 0:
        return f"<p><em>No heatmap for {_esc(title)}</em></p>"
    cw = max(4, min(14, 900 // cols))
    ch = max(4, min(14, 600 // rows))
    pad_l, pad_t, pad_r, pad_b = 90, 28, 30, 18
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
            parts.append(f'<rect x="{pad_l+j*cw}" y="{y}" width="{cw}" height="{ch}" '
                         f'fill="{_heat_color(v, vmax)}" stroke="#fff" stroke-width="0.2">'
                         f'<title>{_esc(s)} / gene{_esc(g)}: {("%.3f"%v) if v is not None else "NA"}</title></rect>')
    # mark filtered gene columns
    for j, g in enumerate(genes):
        if g in filt_genes:
            x = pad_l + j * cw + cw / 2
            parts.append(f'<polygon points="{x-3},{pad_t+rows*ch+2} {x+3},{pad_t+rows*ch+2} {x},{pad_t+rows*ch+8}" fill="red"/>')
    parts.append("</svg>")
    return "\n".join(parts)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strains", required=True, help="protal strain output dir")
    ap.add_argument("--qcmsa", default=None, help="dir with qcmsa outputs (default: --strains)")
    ap.add_argument("--out", required=True, help="report output dir")
    args = ap.parse_args(argv)
    qcmsa_dir = args.qcmsa or args.strains
    os.makedirs(args.out, exist_ok=True)

    species = discover_species(args.strains)
    rows, checks = [], []

    def check(level, sp, msg):
        checks.append((level, sp, msg))

    heat_data = []  # (sp, cell, samples, genes, filt_samples, filt_genes, vmax)

    for sp in species:
        meta = os.path.join(args.strains, sp + ".meta.tsv")
        base_part = os.path.join(args.strains, sp + ".partition.txt")
        pergene = os.path.join(args.strains, sp + ".pergene_filtered.msa.fna")
        base_msa = os.path.join(args.strains, sp + ".msa.fna")
        before_msa = pergene if os.path.exists(pergene) else base_msa

        meta_samples, meta_genes, cell = set(), set(), {}
        vmax = 0.0
        for r in read_tsv(meta):
            meta_samples.add(r["sample"]); meta_genes.add(r["gene_id"])
            try:
                v = float(r.get("multi_rate_vcov2", 0) or 0)
            except ValueError:
                v = 0.0
            cell[(r["sample"], r["gene_id"])] = v
            vmax = max(vmax, v)

        b_seqs, b_samp = count_fasta_seqs(before_msa)
        b_genes = count_partition_genes(base_part)
        summ = load_summary(os.path.join(qcmsa_dir, sp + ".qcmsa_summary.tsv"))
        a_seqs, a_samp = count_fasta_seqs(os.path.join(qcmsa_dir, sp + ".filtered.msa.fna"))
        a_genes = count_partition_genes(os.path.join(qcmsa_dir, sp + ".filtered.partition.txt"))
        snp_agg, snp_n = aggregate_snp_stats(os.path.join(args.strains, sp + ".snp_stats.tsv"))

        row = dict(species=sp, meta_samples=len(meta_samples), meta_genes=len(meta_genes),
                   before_seqs=b_seqs, before_sample_seqs=b_samp, before_genes=b_genes,
                   after_seqs=a_seqs, after_sample_seqs=a_samp, after_genes=a_genes,
                   qcmsa_ran=summ is not None)
        if summ:
            for k in ("sites_in", "sites_kept", "genes_filtered", "samples_filtered", "outlier_cells"):
                row[k] = summ["count"].get(k)
        rows.append(row)

        if b_seqs and len(meta_samples) >= 2 and len(meta_genes) >= 2:
            fs = {s for s, _, _ in summ["samples"]} if summ else set()
            fg = {g for g, _, _ in summ["genes"]} if summ else set()
            heat_data.append((sp, cell, sorted(meta_samples),
                              sorted(meta_genes, key=lambda g: int(g) if g.isdigit() else g),
                              fs, fg, vmax))

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
    write_html(args.out, rows, checks, species, args.strains, qcmsa_dir, heat_data)

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


def write_html(out, rows, checks, species, strains_dir, qcmsa_dir, heat_data):
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

    # ---- heatmaps ----
    heat_svgs = [svg_heatmap(s, g, cell, fs, fg,
                             f"{sp.replace('s__','')}  -  MRate2 (red label/▲ = qcmsa-filtered)",
                             max(0.02, vmax))
                 for (sp, cell, s, g, fs, fg, vmax) in heat_data]

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
        trs.append(
            f"<tr><td>{_esc(r['species'])}</td><td>{r['meta_samples']}</td>"
            f"<td>{r['before_seqs']}</td><td>{r['before_genes']}</td>"
            f"<td>{r.get('after_seqs','—')}</td><td>{r.get('after_genes','—')}</td>"
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

<h2>Samples &amp; genes per MSA (protal default &rarr; after qcmsa)</h2>
<table>
<tr><th>species</th><th>meta samples</th><th>before seqs</th><th>before genes</th>
<th>after seqs</th><th>after genes</th><th>sites in&rarr;kept</th>
<th>genes filt</th><th>samples filt</th></tr>
{''.join(trs)}
</table>
<p class="muted">No protal MSA (too few samples/genes to reconstruct a strain): {_esc(', '.join(no_msa) or '(none)')}</p>

<h2>(a) Which SNPs were filtered out, and why</h2>
<p class="muted">Per species, summed across samples. Variants are filtered by protal's M1 SNP gates:
low cumulative phred-sum (<code>--snp_min_phred_sum</code>) or insufficient supporting reads
(<code>--snp_min_cov</code>). The second panel shows why positions had no callable variant.</p>
<h3>Variant fate (retained vs filtered by M1 gates)</h3>
<div class="grid">{svg_a}</div>
<h3>Why positions had no callable variant</h3>
<div class="grid">{svg_a2}</div>

<h2>(b) Gene &amp; sample filtering (qcmsa, milestone M5)</h2>
<p class="muted">qcmsa removes genes/samples that are multi-allelicity (MRate2) outliers via the
iterative Tukey-IQR rule. Bars show kept vs filtered; the heatmaps below show the underlying
per-cell MRate2 with filtered rows (red labels) and genes (red &#9650;) flagged.</p>
<h3>Gene filtering</h3>
<div class="grid">{svg_b_genes}</div>
<h3>Sample filtering</h3>
<div class="grid">{svg_b_samples}</div>
{''.join('<div class="grid">'+h+'</div>' for h in heat_svgs)}

<h2>Filtering reasons (detail)</h2>
{''.join(reason_blocks) or '<p class="muted">No genes or samples were filtered by qcmsa on this run.</p>'}

<p class="muted">Generated by scripts/strain_test/strain_report.py</p>
</body></html>"""
    with open(os.path.join(out, "report.html"), "w") as fh:
        fh.write(doc)


if __name__ == "__main__":
    sys.exit(main())
