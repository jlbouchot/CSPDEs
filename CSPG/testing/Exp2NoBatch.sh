#!/bin/bash
#title          :Exp2NoBatch.sh
#description    :This script tests the influence of the number of parameters keeping fixed all other parameters.
#author         :Jean-Luc Bouchot
#date           :2026/11/09 [Completely overwriting previous versions]
#updates        :[]
#version        :0.1  
#usage          :bash Exp2NoBatch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 2: h final fixed, number of multi-levels fixed. 
# 1. Call the small driver_MLCSPG_modulated_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_modulated_2D.py script with more options making sure everything works fine.

outdir="results/dgrows"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    d_values=(3 5 7 10 13 16 20 24 28)
    no_compute=False
else
    # Debug execution.
    d_values=(3 5 7)
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

outdir=${outdir}V2
logdir="$outdir/logs"
# Execute the python command.
mkdir -p "$outdir" ### TODO: Review this once agreed upon some naming conventions
mkdir -p "$logdir"
for i in $(seq 0 $((${#d_values[@]} - 1))); do
  # Define output file name based on the current h_0 value.
  output_file="d_${d_values[$i]}"
  # Define a clear log file
  log_file="$logdir/dgrows_${d_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  python driver_MLCSPG_modulated_2D.py --cfg ../data/Exp2_dvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_cosines ${d_values[$i]} --output_file "$output_file" > "$log_file"
done
