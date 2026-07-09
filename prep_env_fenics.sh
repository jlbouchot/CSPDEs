#!/bin/bash
# File prep_env_fenics.sh
# Version 1.0
# 2026/07/06, Jean-Luc Bouchot, jean-luc.bouchot@inria.fr
conda create --name env-fenics36 python=3.6
conda activate env-fenics36
conda install -c conda-forge fenics
# conda install -c conda-forge numpy=1.21.6
conda install -c conda-forge superlu_dist=6.2.0
# Test things
## Dolfin et al
python -c "from dolfin import *; info(NonlinearVariationalSolver.default_parameters(), True)"
## MLCSPG things
conda install cvxpy==1.1
conda install numba



### Not used yet but could be important

#### For MLCSPG
conda install progressbar2
conda install matplotlib

#### For fenics and superlu dist
conda install -c conda-forge boost-cpp=1.72.0
conda install -c conda-forge mpi4py=3.0

