#!/bin/bash
#
#SBATCH --job-name=Exp3MLCSPG
#SBATCH --output=res_Exp3_job_%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --partition=gamma,normal
#SBATCH --error error_Exp3_%A_%a.out
#
#SBATCH --array=0-6

#title          :Exp3Batch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/08/18 [Created 2026/08/18]
#updates        :[]
#version        :0.1   
#usage          :bash Exp3Batch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 3: h final fixes, L fixed, h0 fixed,number of multi-levels increasing. 
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp3"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    J_values=(6 5 4 3 2 1 0)
    L_val=6
    tol_res_val=0.00003125
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
mpirun python driver_MLCSPG_2D.py --cfg ../data/Exp3_Jvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level $L_val --l_start ${J_values[$SLURM_ARRAY_TASK_ID]} --output_file "J_val_${J_values[$SLURM_ARRAY_TASK_ID]}" -e $tol_res_val

