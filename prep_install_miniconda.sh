#!/bin/bash
# File prep_install_miniconda.sh
# Version 1.0
# 2026/07/06, Jean-Luc Bouchot, jean-luc.bouchot@inria.fr
mkdir -p /$SCRATCH/$USER/miniconda3
ln -s /$SCRATCH/$USER/miniconda3 $HOME/.local/miniconda3
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh -O install_miniconda.sh
bash install_miniconda.sh -b -p $HOME/.local/miniconda3
# Activate conda env by default, assuming you are using bash
$HOME/.local/miniconda3/bin/conda init bash
