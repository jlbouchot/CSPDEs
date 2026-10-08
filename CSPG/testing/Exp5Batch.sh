#!/bin/bash
#
#SBATCH --job-name=Exp3MLCSPG
#SBATCH --output=res_Exp3_job_%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --partition=gamma,normal
#SBATCH --error error_Exp3_%A_%a.out
#
#SBATCH --array=0-2

#title          :Exp5Batch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/08/20 [Created 2026/08/20]
#updates        :[]
#version        :0.1   
#usage          :bash Exp5Batch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 5: everything is fixed but the number of cosine dimensions d
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp5"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    # h_0_values=(640 320 160 80 40 20 10)
    d_values=(7 10 13 16 20 24 28 33 38)
else
    # Debug execution.
    d_values=(7 10 13)
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
mpirun python driver_MLCSPG_2D.py --cfg ../data/Exp5_dvaries.ini --experiment_name "$outdir" --no_compute $no_compute --output_file "d_val_${d_values[$SLURM_ARRAY_TASK_ID]}" --nb_cosines ${d_values[$SLURM_ARRAY_TASK_ID]}"
