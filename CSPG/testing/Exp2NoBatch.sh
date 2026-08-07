#!/bin/bash
#title          :Exp2NoBatch.sh
#description    :This script tests the influence of the starting h0 for a given target accuracy and number of approximation levels.
#author         :Jean-Luc Bouchot
#date           :2026/08/04 [Created 2026/07/23]
#updates        :[0.2: 2026/08/04: Updated the script to include a debug mode for testing with smaller values of h0 and fewer levels, updating tol res values.]
#version        :0.2   
#usage          :bash Exp2NoBatch.sh
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


logdir="$outdir/logs"
# Execute the python command.
mkdir -p "$outdir" ### TODO: Review this once agreed upon some naming conventions
mkdir -p "$logdir"
for i in $(seq 0 $((${#h_0_values[@]} - 1))); do
  # Define output file name based on the current h_0 value.
  # output_file="$outdir/test_hFinalFixed_h0_${h_0_values[$i]}"
  output_file="h0_${h_0_values[$i]}"
  # Define a clear log file
  log_file="$logdir/hFinalFixed_h0_${h_0_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  python driver_MLCSPG_2D.py --mesh_x ${h_0_values[$i]} --mesh_y ${h_0_values[$i]} --cfg ../data/Exp2_h0varies.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level ${L_values[$i]} --l_start ${J_values[$i]} --output_file "$output_file" > "$log_file"
done
