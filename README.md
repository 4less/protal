# Protal

Protal is a reference-based taxonomic profiler for bacterial communities and uses paired-end short reads from shotgun metagenomic sequencing as an input. The index is prebuilt and covers the whole taxonomic space from GTDB version r207. The index will soon be made available for download.

# Installation
Protal is in the final steps of development and will then be made available via conda. In the meantime, you can use a local build process via conda as described below.

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
mkdir conda_build
conda build conda_recipe -c conda-forge --output-folder conda_build

# If everything is successful, the local conda package is here
conda_build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## 4. Install in conda 

```{r bash}
# Current directory is your local clone of this repository
conda create -n protal_env conda_build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## Test the installation

```{r bash}
conda activate protal_env
protal
```
