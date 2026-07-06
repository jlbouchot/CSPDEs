#!/bin/bash
# File prep_env_fenics.sh
# Version 1.0
# 2026/07/06, Jean-Luc Bouchot, jean-luc.bouchot@inria.fr
conda create --name env-fenics36 python=3.6
conda activate env-fenics36
conda install -c conda-forge fenics
conda install cvxpy
conda install progressbar2
conda install numba
conda install matplotlib
# conda install -c conda-forge numpy=1.21.6
conda install -c conda-forge superlu_dist=6.2.0
conda install -c conda-forge boost-cpp=1.72.0
conda install -c conda-forge mpi4py=3.0
