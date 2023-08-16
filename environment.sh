#!/usr/bin/bash

CONDA_PREFIX=/usr/users/QIB_fr017/fritsche/miniconda3/envs/gxx
GCC_PREFIX=/usr/users/QIB_fr017/fritsche/miniconda3/envs/gxx/x86_64-conda-linux-gnu/

export COMPILER_PATH=${CONDA_PREFIX}/bin:${GCC_PREFIX}/bin:$COMPILER_PATH
export LIBRARY_PATH=${CONDA_PREFIX}/lib:${GCC_PREFIX}/lib:$LIBRARY_PATH
export C_INCLUDE_PATH=${CONDA_PREFIX}/include:${GCC_PREFIX}/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=${CONDA_PREFIX}/include:${GCC_PREFIX}/include:$CPLUS_INCLUDE_PATH

echo $COMPILER_PATH
echo $LIBRARY_PATH
echo $C_INCLUDE_PATH
echo $CPLUS_INCLUDE_PATH

GOOFY="hello"

export GOOFY