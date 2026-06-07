set shell := ["bash", "-cu"]

# ---- strain test (M1-M5) parameters -- override on the command line, e.g.
#   just strain-test preset=sensitive
# ----------------------------------------------------------------------------
protal      := "build/protal"
strain_db   := env_var_or_default("PROTAL_DB_PATH", "/home/fritscher/data/db/tool/protal/protal-db-r226-0.5.1a")
strain_test2:= "/home/fritscher/non-git/simulate_metagenomes_test/test2"
strain_run  := justfile_directory() / "strain_test_out"
strain_threads := "8"
preset      := "default"

# Delete all cmake-build-
clear:
    rm -rf cmake-build-*

# Full M1-M5 strain test: run protal (default settings, with the qcmsa post-filter)
# on the pre-existing test2 alignments, then build the HTML QC report.
strain-test: strain-protal strain-dbcounts strain-report
    @echo "Report: {{strain_run}}/report/report.html"

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
      printf '#INPUT_DIR\t%s\n'           "{{strain_test2}}/output/reads"; \
      printf '#SAM_OUTPUT_DIR\t%s\n'      "{{strain_test2}}/protal/alignments"; \
      printf '#PROFILE_OUTPUT_DIR\t%s\n'  "{{strain_run}}/profiles"; \
      printf '#STRAIN_OUTPUT_DIR\t%s\n'   "{{strain_run}}/strains"; \
      printf '#MISC_OUTPUT_DIR\t%s\n'     "{{strain_run}}/misc"; \
      grep -vE '^#(INPUT_DIR|OUTPUT_DIR)' "{{strain_test2}}/protal_map.tsv"; \
    {{ '}' }} > {{strain_run}}/strain_test_map.tsv
    PROTAL_DB_PATH="{{strain_db}}" {{protal}} profile \
        --map {{strain_run}}/strain_test_map.tsv \
        -t {{strain_threads}} \
        --run_qcmsa --strain_preset {{preset}} \
        2>&1 | tee {{strain_run}}/protal_run.log | tr '\r' '\n' | grep -vE '^\[=*>* *\] *[0-9]+ %' || true

# Build a per-species ML tree from each strain MSA with IQ-TREE. Override the
# input with strain_msa_suffix=pergene_filtered (protal, more sites) or =filtered
# (qcmsa final, default). Uses the `iqtree` conda env by default; override the
# launcher with strain_iqtree=... . Skips gracefully if IQ-TREE is unavailable.
strain_msa_suffix := "filtered"
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
    built=0
    for msa in "{{strain_run}}/strains/"*."{{strain_msa_suffix}}".msa.fna; do
        sp=$(basename "$msa" ."{{strain_msa_suffix}}".msa.fna)
        nseq=$(grep -c '^>' "$msa")
        if [ "$nseq" -lt 4 ]; then
            echo "[trees] $sp: only $nseq sequences (<4) - skipping"
            continue
        fi
        echo "[trees] $sp ($nseq taxa) -> {{strain_run}}/trees/$sp.treefile"
        # Unpartitioned GTR+G ML tree with 1000 ultrafast bootstraps. IQ-TREE
        # reads the IUPAC ambiguity codes protal writes (milestone M2) natively.
        # A matching .{{strain_msa_suffix}}.partition.txt is available for
        # partitioned analyses (add: -p "$part").
        $iq -s "$msa" -m GTR+G -B 1000 -T AUTO --seqtype DNA \
            --prefix "{{strain_run}}/trees/$sp" -redo \
            > "{{strain_run}}/trees/$sp.iqtree.stdout.log" 2>&1 \
          && built=$((built+1)) \
          || echo "[trees] $sp: IQ-TREE failed (see {{strain_run}}/trees/$sp.iqtree.stdout.log)"
    done
    echo "[trees] built $built tree(s) in {{strain_run}}/trees (suffix=.{{strain_msa_suffix}})"

# (Re)build the self-contained HTML QC report from an existing strain run.
strain-report:
    python3 scripts/strain_test/strain_report.py \
        --strains {{strain_run}}/strains \
        --out {{strain_run}}/report
    @echo "Open {{strain_run}}/report/report.html"

# Remove the strain test outputs.
strain-clean:
    rm -rf {{strain_run}}

# Baseline (no AVX) build
baseline:
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target protal -- -j$(nproc)

# AVX2 build
avx2:
    cmake -S . -B cmake-build-avx2 -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-avx2 --target protal_avx2 -- -j$(nproc)

# Baseline simulate_metagenomes build
simulate:
    cmake -S . -B cmake-build-baseline -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-baseline --target simulate_metagenomes -- -j$(nproc)

# Static baseline build (protal_static target)
static:
    cmake -S . -B cmake-build-static -DCMAKE_BUILD_TYPE=Release
    cmake --build cmake-build-static --target protal_static -- -j$(nproc)

# Build all targets, baseline avx2 static and simulate
build-all: clear baseline avx2 static simulate

# Install protal, protal_avx2, protal_map_utils, protal_launcher and simulate_metagenomes into prefix/bin
install prefix="$HOME/.local":
    mkdir -p {{prefix}}/bin
    cp cmake-build-baseline/protal                  {{prefix}}/bin/protal_baseline
    cp cmake-build-avx2/protal_avx2                 {{prefix}}/bin/protal_avx2
    cp cmake-build-baseline/simulate_metagenomes    {{prefix}}/bin/simulate_metagenomes
    cp scripts/protal_map_utils                     {{prefix}}/bin/protal_map_utils
    cp scripts/protal_launcher                      {{prefix}}/bin/protal
    chmod +x {{prefix}}/bin/protal_baseline {{prefix}}/bin/protal_avx2 {{prefix}}/bin/simulate_metagenomes {{prefix}}/bin/protal_map_utils {{prefix}}/bin/protal
