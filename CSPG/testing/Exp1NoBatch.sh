#!/bin/bash
#title          :Exp1NoBatch.sh
#description    :This script tests the influence of the target accuracy (as illustrated by the number of levels) for a given starting h0 for a given number of approximation levels.
#author         :Jean-Luc Bouchot
#date           :2026/11/04 [Created 2026/08/04]
#updates        :[Checks for the values of the various constants in the process to validate sizes and constants.]
#version        :0.2   
#usage          :bash Exp1NoBatch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 1: h final evolving, number of multi-levels increasing, h0 fixed. 
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/h0_nbDisc_Fixed"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    J_values=(0 1 2 3 4 5 6)
    L_values=(2 3 4 5 6 7 8)
else
    # Debug execution.
    J_values=(0 1 2)
    L_values=(2 3 4)
    outdir="${outdir}_debug"
fi

if [[ -z "$NO_COMPUTE" ]]; then
    # Only check for the values of the various constants in the process to validate sizes and constants.
    no_compute=False
else
    no_compute=True
    outdir="${outdir}_nocompute"
fi

outdir=${outdir}V2
logdir="$outdir/logs"
# Execute the python command.
mkdir -p "$outdir" ### TODO: Review this once agreed upon some naming conventions
mkdir -p "$logdir"
for i in $(seq 0 $((${#L_values[@]} - 1))); do
  # Define output file name based on the current h_0 value.
  output_file="L_val_${L_values[$i]}"
  # Define a clear log file
  log_file="$logdir/h0FixedLfinal${L_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  # python driver_MLCSPG_modulated_2D.py --cfg ../data/Exp1_hfvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level ${L_values[$i]} --l_start ${J_values[$i]} --output_file "$output_file" -e ${tol_res_values[$i]} > "$log_file"
  python driver_MLCSPG_modulated_2D.py --cfg ../data/Exp1_hfvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level ${L_values[$i]} --l_start ${J_values[$i]} --output_file "$output_file" > "$log_file"
done