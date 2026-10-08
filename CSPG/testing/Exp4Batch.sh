#!/bin/bash
#
#SBATCH --job-name=Exp4MLCSPG
#SBATCH --output=res_Exp4_job_%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --partition=gamma,normal
#SBATCH --error error_Exp4_%A_%a.out
#
#SBATCH --array=0-6

#title          :Exp4Batch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/08/18 [Created 2026/08/18]
#updates        :[]
#version        :0.1   
#usage          :bash Exp4Batch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 4: h final fixes, L fixed, h0 fixed,number of multi-levels increasing. 
# This is in fact the same as EXP number 3 but with "easier" values to compute.
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp4"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    # h_0_values=(640 320 160 80 40 20 10)
    J_values=(6 5 4 3 2 1 0)
    L_val=6
    tol_res_val=0.0000078125
else
    # Debug execution.
    J_values=(4 3 2)
    L_val=4
    tol_res_val=0.000125
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
mpirun python driver_MLCSPG_2D.py --cfg ../data/Exp4_Jvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level $L_val --l_start ${J_values[$SLURM_ARRAY_TASK_ID]} --output_file "J_val_${J_values[$SLURM_ARRAY_TASK_ID]}" -e $tol_res_val

