#!/bin/bash
#title          :Exp5NoBatch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/08/18 [Created 2026/08/18]
#updates        :[]
#version        :0.1   
#usage          :bash Exp5NoBatch.sh
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


logdir="$outdir/logs"
# Execute the python command.
mkdir -p "$outdir" ### TODO: Review this once agreed upon some naming conventions
mkdir -p "$logdir"
for i in $(seq 0 $((${#d_values[@]} - 1))); do
  # Define output file name based on the current d value.
  # output_file="$outdir/test_d_${d_values[$i]}"
  output_file="d_val_${d_values[$i]}"
  # Define a clear log file
  log_file="$logdir/d${d_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  python driver_MLCSPG_2D.py --cfg ../data/Exp5_dvaries.ini --experiment_name "$outdir" --no_compute $no_compute --output_file "$output_file" --nb_cosines ${d_values[$i]} > "$log_file"
done
