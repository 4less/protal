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
  --genome-table genomes.tsv \
  --output-dir sims/ \
  --samples 3 \
  --sample-prefix sim \
  --total-read-pairs 100000 \
  --genomes-per-sample 15 \
  --distribution power_law \
  --strains-per-species "Escherichia coli=2,Bacillus subtilis=1"
```
Reads are simulated with `art_illumina`, concatenated per sample into `<sample>_R1.fq` and `<sample>_R2.fq`, and a `manifest.tsv` records the composition.
