#!/bin/bash
# GAMMCOR environment setup — source this file or it's sourced by other scripts
# Usage: source env.sh

export MKL_ROOT=/opt/intel/oneapi/mkl/2025.3
export LD_LIBRARY_PATH=${MKL_ROOT}/lib:${HOME}/.local/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${HOME}/.local/lib/pkgconfig:${PKG_CONFIG_PATH}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
