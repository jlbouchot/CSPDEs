#!/bin/bash
#title          :Exp3NoBatch.sh
#description    :This script tests the influence of the number of discretization levels keeping fixed the starting h0 and the target accuracy.
#author         :Jean-Luc Bouchot
#date           :2026/11/09 [Created 2026/08/04]
#updates        :[0.2: Moving to the new version of the driver_MLCSPG_modulated_2D.py]
#version        :0.2   
#usage          :bash Exp3NoBatch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 3: h final fixes, L fixed, h0 fixed,number of multi-levels increasing. 
# 1. Call the small driver_MLCSPG_modulated_2D.py script with its specific config file and redirecting logs
# 2. Call the small driver_MLCSPG_modulated_2D.py script with more options making sure everything works fine.

outdir="results/Exp3Jvaries"
if [[ -z "$DEBUG_MODE" ]]; then
    # Default execution.
    J_values=(6 5 4 3 2)
    L_val=6
else
    # Debug execution.
    J_values=(4 3 2)
    L_val=4
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
for i in $(seq 0 $((${#J_values[@]} - 1))); do
  # Define output file name based on the current h_0 value.
  output_file="J_val_${J_values[$i]}"
  # Define a clear log file
  log_file="$logdir/h0hfFixedJstart${J_values[$i]}"
  # Execute the python command with the current h_0 value and redirect the output to the corresponding log file.
  python driver_MLCSPG_modulated_2D.py --cfg ../data/Exp3_Jvaries.ini --experiment_name "$outdir" --no_compute $no_compute --nb_level $L_val --l_start ${J_values[$i]} --output_file "$output_file" > "$log_file"
done
