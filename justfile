set shell := ["bash", "-cu"]

# ---- strain test (M1-M5) parameters -- override on the command line, e.g.
#   just strain-test preset=sensitive
# ----------------------------------------------------------------------------
build_dir   := "build"
protal      := build_dir / "protal"
strain_db   := env_var_or_default("PROTAL_DB_PATH", "/home/fritscher/data/db/tool/protal/protal-db-r226-0.5.1a")
strain_input:= "/home/fritscher/non-git/simulate_metagenomes_test/test2"   # simulated dataset (input)
strain_out  := justfile_directory() / "strain_test_out"
strain_variant := "test1"                 # output subfolder: test1 (full filtering) / test2 (raw)
strain_run  := strain_out / strain_variant
strain_threads := "8"
preset      := "default"
# Extra protal flags. test2 (raw) disables ALL protal gene/sample filtering so the
# MSA keeps every observed gene and detected sample (only M1/M2 SNP filtering stays),
# leaving gene/sample filtering entirely to qcmsa and letting users re-filter.
strain_protal_filter_args := ""

# Delete all build trees
clear:
    rm -rf cmake-build-* {{build_dir}}

# Full M1-M5 strain test: run protal (default settings, with the qcmsa post-filter)
# on the pre-existing dataset alignments, then build the HTML QC report.
strain-test: strain-protal strain-dbcounts strain-report
    @echo "Report: {{strain_run}}/report/report.html"

# "Raw" variant -> strain_test_out/test2: protal filters ONLY SNPs (M1/M2); it does
# NOT drop genes (M3), samples (msa_min_samples / per-seq hcov) or mask by
# multi-allelicity (M4). qcmsa then does all gene/sample filtering. The unfiltered
# <species>.raw.msa.fna can be re-filtered with other thresholds.
strain-test-raw:
    just strain_variant=test2 \
         strain_protal_filter_args="--gene_min_hcov_frac 0 --gene_min_mean_depth 0 --msa_min_samples 0 --msa_min_hcov 0" \
         strain-test

# Count marker genes per species in the DB genome (the true gene denominator,
# revealing how many markers were lost to abundance before M3). Cached as a TSV.
strain-dbcounts:
    PROTAL_DB_PATH="{{strain_db}}" python3 scripts/strain_test/db_gene_counts.py \
        --db "{{strain_db}}" \
        --strains {{strain_run}}/strains \
        --out {{strain_run}}/db_gene_counts.tsv

# Build a map that reuses the existing alignments + reads and writes strain
# results into the repo-local run dir, then run protal with --run_qcmsa (M5).
strain-protal:
    mkdir -p {{strain_run}}/strains {{strain_run}}/misc {{strain_run}}/profiles
    # Assemble the map: base OUTPUT_DIR + reuse existing SAMs/reads, local strain output.
    {{ '{' }} \
      printf '#OUTPUT_DIR\t%s\n'          "{{strain_run}}"; \
      printf '#INPUT_DIR\t%s\n'           "{{strain_input}}/output/reads"; \
      printf '#SAM_OUTPUT_DIR\t%s\n'      "{{strain_input}}/protal/alignments"; \
      printf '#PROFILE_OUTPUT_DIR\t%s\n'  "{{strain_run}}/profiles"; \
      printf '#STRAIN_OUTPUT_DIR\t%s\n'   "{{strain_run}}/strains"; \
      printf '#MISC_OUTPUT_DIR\t%s\n'     "{{strain_run}}/misc"; \
      grep -vE '^#(INPUT_DIR|OUTPUT_DIR)' "{{strain_input}}/protal_map.tsv"; \
    {{ '}' }} > {{strain_run}}/strain_test_map.tsv
    # Write the full log to a file (avoids a pipeline whose trailing grep can fail
    # the recipe under `set -o pipefail`, which lmod's BASH_ENV enables), then show
    # a progress-bar-stripped tail.
    PROTAL_DB_PATH="{{strain_db}}" {{protal}} profile \
        --map {{strain_run}}/strain_test_map.tsv \
        -t {{strain_threads}} \
        --run_qcmsa --strain_preset {{preset}} {{strain_protal_filter_args}} \
        > {{strain_run}}/protal_run.log 2>&1 || true
    -tr '\r' '\n' < {{strain_run}}/protal_run.log | grep -vE '^\[=*>* *\] *[0-9]+ %' | tail -40

# Build a per-species ML tree from each strain MSA with IQ-TREE. strain_tree_input
# selects the MSA: "filtered" = qcmsa output <sp>.msa.fna (default); "raw" = protal
# native <sp>.raw.msa.fna. Uses the `iqtree` conda env by default; override the
# launcher with strain_iqtree=... . Skips gracefully if IQ-TREE is unavailable.
strain_tree_input := "filtered"
strain_iqtree := ""
strain-trees:
    #!/usr/bin/env bash
    set -uo pipefail
    # Resolve an IQ-TREE launcher: explicit override > PATH > `iqtree` conda env.
    iq="{{strain_iqtree}}"
    if [ -z "$iq" ]; then
        iq=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
    fi
    if [ -z "$iq" ] && command -v conda >/dev/null 2>&1; then
        for b in iqtree3 iqtree2 iqtree; do
            if conda run -n iqtree "$b" --version >/dev/null 2>&1; then
                iq="conda run -n iqtree $b"; break
            fi
        done
    fi
    if [ -z "$iq" ]; then
        echo "[trees] IQ-TREE not found (PATH or 'iqtree' conda env); skipping."
        exit 0
    fi
    echo "[trees] using: $iq"
    mkdir -p "{{strain_run}}/trees"
    shopt -s nullglob
    case "{{strain_tree_input}}" in raw) ext=".raw.msa.fna";; *) ext=".msa.fna";; esac
    built=0
    for meta in "{{strain_run}}/strains/"*.meta.tsv; do
        sp=$(basename "$meta" .meta.tsv)
        msa="{{strain_run}}/strains/$sp$ext"
        [ -f "$msa" ] || continue
        nseq=$(grep -c '^>' "$msa")
        if [ "$nseq" -lt 4 ]; then
            echo "[trees] $sp: only $nseq sequences (<4) - skipping"
            continue
        fi
        echo "[trees] $sp ($nseq taxa) -> {{strain_run}}/trees/$sp.treefile"
        # Unpartitioned GTR+G ML tree with 1000 ultrafast bootstraps. IQ-TREE
        # reads the IUPAC ambiguity codes protal writes (milestone M2) natively.
        $iq -s "$msa" -m GTR+G -B 1000 -T AUTO --seqtype DNA \
            --prefix "{{strain_run}}/trees/$sp" -redo \
            > "{{strain_run}}/trees/$sp.iqtree.stdout.log" 2>&1 \
          && built=$((built+1)) \
          || echo "[trees] $sp: IQ-TREE failed (see {{strain_run}}/trees/$sp.iqtree.stdout.log)"
    done
    echo "[trees] built $built tree(s) in {{strain_run}}/trees (input={{strain_tree_input}})"

# Re-filter a raw run's MSAs with qcmsa, applying M3-equivalent coverage gating
# (from the meta hcov/depth columns) PLUS the usual MRate2 + site cleanup. Lets
# you re-filter strain_test_out/test2 (the raw run) with any thresholds without
# re-running protal. Outputs into <run>/refiltered/.
refilter_hcov        := "0.3"
refilter_depth       := "1"
refilter_min_samples := "3"
refilter_sample_abs  := "0"    # >0: remove a sample multi-allelic in >= N genes (catches conspecific/mixed strains)
refilter_gene_abs    := "0"    # >0: remove a gene multi-allelic in >= N samples
strain-refilter:
    #!/usr/bin/env bash
    set -uo pipefail
    mkdir -p "{{strain_run}}/refiltered"
    shopt -s nullglob
    n=0
    for msa in "{{strain_run}}/strains/"*.raw.msa.fna; do
        sp=$(basename "$msa" .raw.msa.fna)
        part="{{strain_run}}/strains/$sp.raw.partition.txt"
        meta="{{strain_run}}/strains/$sp.meta.tsv"
        [ -f "$part" ] && [ -f "$meta" ] || continue
        python3 scripts/qcmsa.py "$msa" "$part" "$meta" \
            --prefix "{{strain_run}}/refiltered/$sp" --preset {{preset}} \
            --gene-min-hcov {{refilter_hcov}} \
            --gene-min-mean-depth {{refilter_depth}} \
            --gene-min-samples {{refilter_min_samples}} \
            --sample-abs-min-bad {{refilter_sample_abs}} \
            --gene-abs-min-bad {{refilter_gene_abs}} \
            > "{{strain_run}}/refiltered/$sp.qcmsa.log" 2>&1 \
          && { echo "[refilter] $sp"; n=$((n+1)); } \
          || echo "[refilter] $sp: failed (see {{strain_run}}/refiltered/$sp.qcmsa.log)"
    done
    echo "[refilter] re-filtered $n species (hcov>={{refilter_hcov}}, depth>={{refilter_depth}}, >{{refilter_min_samples}} samples) -> {{strain_run}}/refiltered"

# (Re)build the self-contained HTML QC report from an existing strain run.
strain-report:
    python3 scripts/strain_test/strain_report.py \
        --strains {{strain_run}}/strains \
        --out {{strain_run}}/report
    @echo "Open {{strain_run}}/report/report.html"

# Remove the strain test outputs.
strain-clean:
    rm -rf {{strain_run}}

# The ISA flags are per-target (isa_baseline / isa_avx2 in CMakeLists.txt), so
# baseline, avx2 and static binaries all come out of one tree -- no need for
# separate cmake-build-* dirs.
# Configure the build tree
configure:
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE=Release

# Baseline (no AVX) build
baseline: configure
    cmake --build {{build_dir}} --target protal -- -j$(nproc)

# AVX2 build
avx2: configure
    cmake --build {{build_dir}} --target protal_avx2 -- -j$(nproc)

# simulate_metagenomes build
simulate: configure
    cmake --build {{build_dir}} --target simulate_metagenomes -- -j$(nproc)

# Static baseline build (protal_static target)
static: configure
    cmake --build {{build_dir}} --target protal_static -- -j$(nproc)

# Build all binaries that `just install` ships
build-all: configure
    cmake --build {{build_dir}} --target protal protal_avx2 simulate_metagenomes -- -j$(nproc)

# Always rebuilds first so the installed binaries match the working tree (an
# out-of-date build dir used to be installed silently).
# Install protal, protal_avx2, protal_map_utils, protal_launcher, qcmsa and simulate_metagenomes into prefix/bin
install prefix="$HOME/.local": build-all
    mkdir -p {{prefix}}/bin
    cp {{build_dir}}/protal                         {{prefix}}/bin/protal_baseline
    cp {{build_dir}}/protal_avx2                    {{prefix}}/bin/protal_avx2
    cp {{build_dir}}/simulate_metagenomes           {{prefix}}/bin/simulate_metagenomes
    cp scripts/protal_map_utils                     {{prefix}}/bin/protal_map_utils
    cp scripts/protal_launcher                      {{prefix}}/bin/protal
    cp scripts/qcmsa.py                             {{prefix}}/bin/qcmsa
    chmod +x {{prefix}}/bin/protal_baseline {{prefix}}/bin/protal_avx2 {{prefix}}/bin/simulate_metagenomes {{prefix}}/bin/protal_map_utils {{prefix}}/bin/protal {{prefix}}/bin/qcmsa
    @echo "Installed to {{prefix}}/bin: $({{prefix}}/bin/protal --version)"
