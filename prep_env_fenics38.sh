#!/bin/bash
# File prep_env_fenics38.sh
# Version 1.0
# 2026/10/08, Jean-Luc Bouchot, jean-luc.bouchot@inria.fr
# Creates the Python 3.8 / FEniCS 2019.1.0 environment for MLCSPG (conda-forge only).
# Usage: bash prep_env_fenics38.sh [--exact]
#   --exact : rebuild the exact verified environment from env-fenics38.lock.txt (linux-64 only)
set -e
ENV_NAME=env-fenics38
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Environment ${ENV_NAME} already exists. Remove it first: conda env remove -n ${ENV_NAME}"
    exit 1
fi

# Mixing conda-forge and defaults breaks legacy FEniCS: force strict channel priority
if [[ "$1" == "--exact" ]]; then
    conda create -y -n ${ENV_NAME} --file "${HERE}/env-fenics38.lock.txt"
else
    CONDA_CHANNEL_PRIORITY=strict conda env create -y -f "${HERE}/env-fenics38.yml"
fi

# Test things: imports, FEniCS solves with the various linear solvers, cvxpy, MPI
conda run --no-capture-output -n ${ENV_NAME} python "${HERE}/check_env_fenics38.py"
conda run --no-capture-output -n ${ENV_NAME} mpirun -n 2 python -c "from dolfin import *; m = UnitSquareMesh(64, 64); print('MPI rank', MPI.rank(m.mpi_comm()), 'cells', m.num_cells())"
echo "Done. Activate with: conda activate ${ENV_NAME}"
