#!/bin/bash
#
#SBATCH --job-name=Exp2MLCSPG
#SBATCH --output=res_Exp2_job_%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --partition=gamma,normal
#SBATCH --error error_Exp2_%A_%a.out
#
#SBATCH --array=0-6

#title          :Exp2Batch.sh
#description    :This script tests the influence of the starting h0 for a given target accuracy and number of approximation levels.
#author         :Jean-Luc Bouchot
#date           :2026/08/18 [Created 2026/08/18]
#updates        :[]
#version        :0.1
#usage          :bash Exp2Batch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 2: h final fixed, number of multi-levels fixed. 
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp2"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    h_0_values=(640 320 160 80 40 20 10)
    J_values=(0 1 2 3 4 5 6)
    L_values=(2 3 4 5 6 7 8)
    no_compute=False
else
    # Debug execution.
    h_0_values=(40 20 10)
    J_values=(0 1 2)
    L_values=(2 3 4)
    no_compute=True
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
mpirun python driver_MLCSPG_2D.py --cfg ../data/Exp2_h0varies.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level ${L_values[$SLURM_ARRAY_TASK_ID]} --l_start ${J_values[$SLURM_ARRAY_TASK_ID]} --mesh_x ${h_0_values[$SLURM_ARRAY_TASK_ID]} --mesh_y ${h_0_values[$SLURM_ARRAY_TASK_ID]} --output_file "h0_${h_0_values[$SLURM_ARRAY_TASK_ID]}"
