# qcmsa — strain MSA post-filtering (how to operate it)

`qcmsa.py` post-filters protal's per-species strain MSAs using the per-`(sample,gene)`
metrics protal writes to `.meta.tsv`. It is the **default** post-filter (protal runs
it automatically) but it is also a standalone tool, so you can **re-filter with any
parameters without re-running protal**.

## Inputs and outputs

For one species, qcmsa reads three files protal produces in the strain output dir:

| input | what it is |
|---|---|
| `<species>.pergene_filtered.msa.fna` | the MSA to filter (use `.msa.fna` if you prefer the unmasked base MSA) |
| `<species>.partition.txt` | gene → column ranges (RAxML style) |
| `<species>.meta.tsv` | per-(sample,gene) coverage + multi-allelicity metrics (the filter's evidence) |

It writes (prefix defaults to the MSA path minus `.fna`):

| output | what it is |
|---|---|
| `<prefix>.filtered.msa.fna` | the filtered MSA |
| `<prefix>.filtered.partition.txt` | partition with recomputed coordinates |
| `<prefix>.qcmsa_summary.tsv` | machine-readable decision log (what was filtered and why) |
| `<prefix>.qc.png` | optional MRate2 heatmap (`--plot`, needs matplotlib) |

## What it filters (four stages, in order)

1. **Coverage gate** (M3-equivalent, **off by default**) — drops genes / gap-fills
   cells that are too sparsely covered, read from the meta `hcov` and
   `mean_vcov_nonzero` columns.
2. **Multi-allelicity (MRate2) filter** — removes genes/samples that are
   multi-allelicity outliers via an iterative Tukey-IQR rule, and masks individual
   outlier cells.
3. **Site cleanup** — drops constant / low-parsimony (uninformative) sites.
4. **reapply-hcov** — drops whole sequences below a valid-base floor.

## The parameters you'll actually change

### Multi-allelicity (contamination / mixed strains)
| flag | default | effect |
|---|---|---|
| `--preset strict\|default\|sensitive` | default | bundle for the two below: strict=1.0/1, default=1.5/2, sensitive=2.0/3 |
| `--iqr-mult FLOAT` | 1.5 | Tukey fence multiplier — **lower = more aggressive** (removes more) |
| `--min-bad INT` | 2 | a gene/sample needs this many bad peers before removal — **higher = less aggressive** |
| `--max-mrate2 FLOAT` | data fence | hard per-cell MRate2 cap for masking outlier cells |
| `--no-mask-cell-outliers` | (masking on) | don't mask individual outlier cells |

### Coverage (M3-equivalent — turn ON to filter low-coverage genes here)
| flag | default | effect |
|---|---|---|
| `--gene-min-hcov FLOAT` | 0 (off) | min fraction of a gene covered for a cell to pass (≈ protal `--gene_min_hcov_frac`) |
| `--gene-min-mean-depth FLOAT` | 0 (off) | min mean depth over covered positions (≈ protal `--gene_min_mean_depth`) |
| `--gene-min-samples INT` | 0 (off) | drop a gene unless **more than** this many samples pass coverage (≈ protal `--msa_min_samples`) |

> To reproduce protal's (new, permissive) M3 inside qcmsa: `--gene-min-hcov 0.3
> --gene-min-mean-depth 1 --gene-min-samples 3`. The old aggressive M3 was
> `0.5 / 3 / 3`.

### Sites and sequences
| flag | default | effect |
|---|---|---|
| `--min-parsimony-samples INT` | 2 | drop sites where fewer than N samples differ from the majority (subsumes constant-site removal) — set `0` to keep all sites |
| `--keep-constant` | (remove on) | keep constant sites |
| `--reapply-hcov INT` | 0 (off) | drop sequences with fewer than N valid bases |

## How to change the parameters — three ways

### 1. Through protal (the easy knob)
protal runs qcmsa automatically. The one parameter it forwards is the preset:

```bash
protal profile --map map.tsv --strain_preset sensitive   # strict | default | sensitive
protal profile --map map.tsv --no_qcmsa                   # skip qcmsa entirely
```

For anything beyond the preset, re-filter the output (below) — you do **not** need
to re-run protal.

### 2. Re-filter an existing run (recommended for tuning)
Run qcmsa directly on the MSA/partition/meta protal already wrote. Nothing is
re-aligned; this takes seconds.

```bash
# stricter multi-allelicity removal + reproduce permissive M3 coverage gating
python3 scripts/qcmsa.py \
    out/strains/s__Bacteroides_ovatus.pergene_filtered.msa.fna \
    out/strains/s__Bacteroides_ovatus.partition.txt \
    out/strains/s__Bacteroides_ovatus.meta.tsv \
    --prefix out/refiltered/s__Bacteroides_ovatus \
    --preset sensitive \
    --gene-min-hcov 0.3 --gene-min-mean-depth 1 --gene-min-samples 3 \
    --min-parsimony-samples 2
```

Inspect `…/s__Bacteroides_ovatus.qcmsa_summary.tsv` to see exactly what was removed
(`genes_filtered_coverage`, `genes_filtered_mrate2`, `samples_filtered`,
`coverage_gap_filled_cells`, `sites_in`→`sites_kept`, and per-gene/sample reasons).

### 3. Batch re-filter via `just` (all species at once)
```bash
just strain_variant=test2 strain-refilter \
     refilter_hcov=0.3 refilter_depth=1 refilter_min_samples=3 preset=sensitive
```
Outputs land in `strain_test_out/<variant>/refiltered/`.

## Picking values — quick guidance

- **Too much gets removed / strains collapse?** Loosen: `--preset sensitive`
  (= higher `--iqr-mult` 2.0 + higher `--min-bad` 3), lower or disable the coverage
  gate, and set `--min-parsimony-samples 0` to keep invariant sites. (Mnemonic:
  higher `--iqr-mult` and higher `--min-bad` = looser; lower = stricter.)
- **Contamination / mixed strains slipping through?** Tighten: `--preset strict`
  (or `--iqr-mult 1.0 --min-bad 1`).
- **Want the rawest possible MSA to feed your own pipeline?** Run protal with
  `--no_qcmsa` (and permissive/zero coverage flags) and filter downstream yourself —
  the `.msa.fna` + `.meta.tsv` carry everything you need.

## Note on where filtering lives

Coverage filtering is **identical** whether protal does it (M3) or qcmsa does it
(`--gene-min-*`) — protal exports the exact stats (`hcov`, `mean_vcov_nonzero`) the
gate uses, so the two produce byte-identical MSAs. Doing it in qcmsa just makes it
**reversible** (re-tune without re-running protal). The **SNP filters (M1)** cannot
move to qcmsa — they need per-read base-quality/strand data that is not in the meta —
so those always stay in protal (`--snp_*`).
