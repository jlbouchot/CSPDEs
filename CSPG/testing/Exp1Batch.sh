#!/bin/bash
#
#SBATCH --job-name=Exp1MLCSPG
#SBATCH --output=res_Exp1_job_%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --partition=gamma,normal
#SBATCH --error error_Exp1_%A_%a.out
#
#SBATCH --array=0-6

#title          :Exp1Batch.sh
#description    :This script tests the influence of the target accuracy (as illustrated by the number of levels) for a given starting h0 for a given number of approximation levels.
#author         :Jean-Luc Bouchot
#date           :2026/08/17 [Created 2026/08/17]
#updates        :[]
#version        :0.1   
#usage          :bash Exp1Batch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 1: h final evolving, number of multi-levels increasing, h0 fixed. 
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp1"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    J_values=(0 1 2 3 4 5 6)
    L_values=(2 3 4 5 6 7 8)
    tol_res_values=(0.0005 0.00025 0.000125 0.0000625 0.00003125 0.000015625 0.0000078125)
else
    # Debug execution.
    J_values=(0 1 2)
    L_values=(2 3 4)
    tol_res_values=(0.0005 0.00025 0.000125)
    outdir="${outdir}_debug"
fi

if [[ -z "$NO_COMPUTE" ]]; then
    # Only check for the values of the various constants in the process to validate sizes and constants.
    no_compute=False
else
    no_compute=True
    outdir="${outdir}_nocompute"
fi

# Execute the python command.
mkdir -p "$outdir" 
mpirun python driver_MLCSPG_2D.py --cfg ../data/Exp1_hfvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level ${L_values[$SLURM_ARRAY_TASK_ID]} --l_start ${J_values[$SLURM_ARRAY_TASK_ID]} --output_file "L_val_${L_values[$SLURM_ARRAY_TASK_ID]}" -e ${tol_res_values[$SLURM_ARRAY_TASK_ID]}
