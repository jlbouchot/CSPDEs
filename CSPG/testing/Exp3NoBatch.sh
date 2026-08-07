#!/bin/bash
#title          :Exp3NoBatch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/08/04 [Created 2026/08/04]
#updates        :[]
#version        :0.1   
#usage          :bash Exp3NoBatch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 3: h final evolving, number of multi-levels increasing, h0 fixed. 
# 1. Call the small driver_MLCSPG_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_2D.py script with more options making sure everything works fine.

outdir="results/Exp3"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    # h_0_values=(640 320 160 80 40 20 10)
    J_values=(6 5 4 3 2 1 0)
    L_val=6
    tol_res_val=0.0000078125
else
    # Debug execution.
    J_values=(3 2 1)
    L_val=3
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


logdir="$outdir/logs"
# Execute the python command.
mkdir -p "$outdir" ### TODO: Review this once agreed upon some naming conventions
mkdir -p "$logdir"
for i in $(seq 0 $((${#J_values[@]} - 1))); do
  # Define output file name based on the current h_0 value.
  # output_file="$outdir/test_hFinalFixed_h0_${h_0_values[$i]}"
  output_file="J_val_${J_values[$i]}"
  # Define a clear log file
  log_file="$logdir/h0hfFixedJstart${J_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  python driver_MLCSPG_2D.py --cfg ../data/Exp3_Jvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level $L_val --l_start ${J_values[$i]} --output_file "$output_file" -e $tol_res_val > "$log_file"
done
