# Protal

Protal is a reference-based taxonomic profiler for bacterial communities and uses paired-end short reads from shotgun metagenomic sequencing as an input. The index is prebuilt and covers the whole taxonomic space from GTDB version r214. The index is available for download under https://protal.earlham.ac.uk/main.php?site=downloads

# Installation
Protal is in the final steps of development and is also available via conda. In the meantime, you can use a local build process via conda as described below.

## Requirements?
- git
- conda
- A linux machine (no support for mac or windows)

## Steps?
1. Install conda-build
2. clone git repository
3. build protal locally with conda-build
4. install in conda environment from local build

## 1. Install conda-build
This is needed to build a conda project from local files.
```{r bash}
conda install conda-build
```

Alternatively, if you are using micromamba or mamba, you can also install conda-build with
```{r bash}
micromamba install conda-build
# or
mamba install conda-build
```

## 2. Clone this repository
Clone this repository.
```{r bash}
git clone git@github.com:4less/protal.git
```

## 3. build protal locally with conda-build
Compiles protal from the source files with instructions supplied in conda-recipe/meta.yml and conda-recipe/build.sh.
```{r bash}
cd protal
mkdir conda-build
conda build conda-recipe -c conda-forge --output-folder conda-build

# If everything is successful, the local conda package is here
conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## 4. Install in conda 

```{r bash}
# Current directory is your local clone of this repository
conda create -n protal_env conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
#or
micromamba create -n protal_env conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## Test the installation

```{r bash}
conda activate protal_env
protal
```

## Metagenome simulation (C++)

Build the simulator helper binary:
```bash
cmake -S . -B cmake-build-release
cmake --build cmake-build-release --target simulate_metagenomes
```

Input TSV format (three columns): genome name, GTDB taxonomy string, path to genome FASTA (supports .gz). Example run:
```bash
./cmake-build-release/simulate_metagenomes \
  --genome_table genomes.tsv \
  --output_dir sims/ \
  --samples 3 \
  --sample_prefix sim \
  --total_read_pairs 100000 \
  --species_per_sample 15 \
  --distribution power_law \
  --strains_per_species "0.4,0.2"
```
Reads are simulated with `art_illumina`, concatenated per sample into `<sample>_R1.fq` and `<sample>_R2.fq`, and a `manifest.tsv` records the composition.

### Reproducing a simulated dataset

Every run writes its full provenance next to the reads:

- `manifest.tsv` — one row per (sample, genome), including the FASTA it came from and
  the `art_seed` ART used for it. `manifests/<sample>.tsv` holds the same rows split
  per sample.
- `run_params.tsv` — the command line, the RNG seed (resolved and recorded even when
  `--seed` was not given), and the ART settings.

A manifest is a complete, self-contained description of a dataset, so replaying one
does not depend on reproducing the community-design RNG:

```bash
./cmake-build-release/simulate_metagenomes \
  --from_manifest sims/manifest.tsv \
  --output_dir sims_replay/
```

The replay reproduces the reads byte for byte. All sampling options
(`--distribution`, `--species_per_sample`, `--seed`, …) are ignored; only the ART
settings still apply, and `--read_length` is checked against the depths the manifest
implies. A per-sample manifest replays just that sample.

Manifests written before the `fasta_path` and `art_seed` columns existed still replay:
pass `--genome_table` so the genomes can be resolved by name, and the composition and
per-genome depth are reproduced exactly while the reads themselves are fresh
realizations. The replay's own manifest carries seeds, so it is exactly reproducible
from then on.
