#!/bin/bash

# Test file to verify that all args are passed downstream and that MLCSPG works accordingly.
# 1. Call the small driver_MLCSPG_1D_runner.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_1D_runner.py script with more options making sure everything works fine.

# Execute the python command.
mkdir -p results/testMLCSPG
python driver_MLCSPG_1D_runner.py --cfg ../data/simple_1d_cfg.ini --output_file wiht > results/testMLCSPG/test_wiht_1D.log
# python driver_MLCSPG_1D_runner.py --cfg ../data/simple_1d_cfg.ini --recovery_algo bpdn --output_file bpdn > results/testMLCSPG/test_bpdn_1D.log # This is already too much for my computer
python driver_MLCSPG_1D_runner.py --cfg ../data/simple_1d_cfg.ini --recovery_algo whtp --output_file whtp > results/testMLCSPG/test_whtp_1D.log
