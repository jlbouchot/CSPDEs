#!/bin/bash
#title          :Exp4NoBatch.sh
#description    :This script tests the influence of the base mesh h0 = 1/mesh_x at a FIXED target accuracy h_L and a fixed number of discretization levels.
#author         :Jean-Luc Bouchot
#date           :2026/11/14 [Created 2026/11/14]
#updates        :[]
#version        :0.1
#usage          :bash Exp4NoBatch.sh
#options        :Pass the env variables DEBUG_MODE for running a small example and/or NO_COMPUTE for only checking the values of the various constants in the process to validate sizes and constants.
#notes          :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Test file to call all tests for the MLCSPG paper, experiment number 4:
# h0 varies, h_L fixed, number of levels fixed.
#
# WHAT THIS ISOLATES
#   The grids are mesh_x * 2^l for l = J..L and the sparsity schedule
#   (CSPDE_v2.level_plan) is a function of (J, L) ALONE.  mesh_x and the level
#   index are therefore two handles on the same mesh, and holding
#       mesh_x * 2^L = 640        (same target mesh on every row)
#       L - J       = 2           (same number of levels on every row)
#   pins the actual grids to 160^2, 320^2, 640^2 for the WHOLE sweep while the
#   schedule slides.  So this experiment changes nothing a PDE solve can see --
#   it moves only the sample budget the theory asks for, which is exactly what
#   makes it a calibration test rather than a convergence study.
#
#   The tolerance is invariant too, which is the point.  CSPDE_v2 derives it
#   from utils.default_tolerance_sparse_recovery(L, L-J+1, mesh_x, t_0, t_prime)
#     = 0.125 / ( mesh_x^2 * 4^L * (2^(t_0+t')-1) * (L-J+1) )
#     = 0.125 / ( (mesh_x * 2^L)^2 * 3 * 3 )
#   and (mesh_x * 2^L) and (L-J+1) are both fixed, so every row recovers to
#   tol_res = 3.39e-08.  Nothing but s_J differs across the sweep.
#
# THE SWEEP  (d = 20; s from level_plan, m = ceil(2 log(s)^2 s log N))
#
#   mesh_x   h0     J  L   grids             plan (level,s)   m per level            PDE
#        5   1/5    5  7   160, 320, 640     264, 26, 13      268746, 5377, 1302     4.2 h
#       10   1/10   4  6   160, 320, 640     133, 26, 13       93177, 5377, 1302     2.0 h
#       20   1/20   3  5   160, 320, 640      68, 26, 13       30950, 5377, 1302     1.3 h
#       40   1/40   2  4   160, 320, 640      51, 26, 13       18953, 5377, 1302     1.1 h
#       80   1/80   1  3   160, 320, 640      51, 26, 13       18953, 5377, 1302     1.1 h
#
#   s_J = max( s_L * 2^((L-J) rate p/(1-p)) , ceil(const_sj^(p0/(1-p0)) 2^(L p0 rate/(1-p0))) )
#   The first term is 51 on every row (it depends on L-J only, which is fixed);
#   the second is the L-dependent const_sj floor: 18, 34, 68, 133, 264 for
#   L = 3..7.  So the floor is what the sweep actually varies, and it takes over
#   at L >= 5.  Hence:
#     * mesh_x = 40 and 80 have IDENTICAL plans on IDENTICAL grids.  They are
#       kept on purpose as a NULL TEST: they must produce the same error and the
#       same cost, and if they do not, something outside level_plan depends on
#       the level labelling.  mesh_x = 160 (J=0, L=2) would be a third copy and
#       is left out.
#     * mesh_x = 5 is a 5x5 nominal base mesh, which is never meshed: only
#       levels J..L = 5..7 are, i.e. 160^2 upwards.  It is a label, not a grid.
#
#   Reading the result: cost rises 3.8x from mesh_x=80 to mesh_x=5 while the
#   grids, the tolerance and the detail levels are untouched.  If the error is
#   flat across that range, the const_sj floor over-provisions the first level
#   and const_sj can be lowered; if the error improves, the floor is load-
#   bearing and Exp1's L-dependence is partly a sparsity effect, not a mesh one.
#
# GROUND TRUTH: nothing to generate, and nothing to wait for.  d, power,
#   weight_cosine, fluctuation_importance, abar, the RHS, the QoI, the GT mesh
#   and the element degree are all shared with Exp1/Exp2/Exp3, and mesh_x is NOT
#   part of the identity (utilities/gt_meta.py), so all five rows score against
#   one reference:
#     qoi = box (this .ini)  GTresults_d20Hrich640-1280-2560_1182f955
#                            -- ALREADY GENERATED and verified
#     qoi = avg (--qoi avg)  GTresults_d20Hrich640-1280-2560_9154097c
#                            -- pending launchGT.sh --qoi avg
#   Exp4 is therefore the one experiment whose figures can be produced as soon
#   as the runs land.  Plot them with utilities/graphWithRespectToH0.py.
#
# 1. Call the driver_MLCSPG_modulated_2D.py script with its specific config file and redirecting logs
# 2. Call it once per (mesh_x, J, L) row of the sweep

set -u

outdir="results/Exp4meshxvaries"

## Parallel arrays, one entry per row. They are checked below.
if [[ -z "${DEBUG_MODE:-}" ]]; then
    # Default execution.
    MX_values=(5 10 20 40 80)
    J_values=(5 4 3 2 1)
    L_values=(7 6 5 4 3)
else
    ## Debug execution: a 160^2 target instead of 640^2, ~2 min per row.
    ## At L = 4, 3, 2 the const_sj floor is 34, 18, 9, all below the 51 that
    ## L-J=2 forces, so ALL THREE ROWS GET THE SAME PLAN (51, 26, 13) on the
    ## same grids (40, 80, 160) with the same tolerance.  That makes the debug
    ## sweep a pure NULL TEST: it must print identical s, N, m and epsilon for
    ## every row, with only the level LABELS shifted (2,3,4 / 1,2,3 / 0,1,2).
    ## Measured: epsilon = 7.469069826324926e-05 on all three, which is also
    ## bit-identical to Exp1's L=4 log -- same target mesh, same level count.
    ## It exercises the mesh_x override, level_plan and the tolerance
    ## derivation without touching a PDE solve.  Run it with NO_COMPUTE=1.
    MX_values=(10 20 40)
    J_values=(2 1 0)
    L_values=(4 3 2)
    outdir="${outdir}_debug"
fi

if [[ -z "${NO_COMPUTE:-}" ]]; then
    # Only check for the values of the various constants in the process to validate sizes and constants.
    no_compute=False
else
    no_compute=True
    outdir="${outdir}_nocompute"
fi

## ------------------------------------------------------------------ guards
## The invariants the experiment rests on. A hand-edited array that breaks one
## of them would still run and still produce plausible numbers, so they are
## checked rather than commented.
TARGET=$(( ${MX_values[0]} * (1 << ${L_values[0]}) ))
NLEVELS=$(( ${L_values[0]} - ${J_values[0]} + 1 ))
for i in "${!MX_values[@]}"; do
    mx=${MX_values[$i]}; J=${J_values[$i]}; L=${L_values[$i]}
    (( mx * (1 << L) == TARGET )) || {
        echo "FATAL: row $i (mesh_x=$mx, L=$L) targets $(( mx * (1 << L) ))^2, not ${TARGET}^2."
        echo "       Every row must share the target mesh or the rows are not comparable."; exit 1; }
    (( L - J + 1 == NLEVELS )) || {
        echo "FATAL: row $i (J=$J, L=$L) has $(( L - J + 1 )) levels, not $NLEVELS."
        echo "       A differing level count changes tol_res and the detail sparsities too."; exit 1; }
    (( J >= 0 )) || { echo "FATAL: row $i has l_start = $J < 0."; exit 1; }
done
## The sampling cache is <experiment_name>/sampling_points_Lmax<L>_d<d>_l<l>_s_<s>.npy
## and does NOT record the mesh. Two rows with the same L would therefore share
## a cache entry while sitting on different grids -- a silent wrong-sample bug.
## Distinct L per row is what rules that out, so it is checked.
if [[ $(printf '%s\n' "${L_values[@]}" | sort | uniq -d | wc -l) -ne 0 ]]; then
    echo "FATAL: repeated nb_level in L_values: $(printf '%s ' "${L_values[@]}")"
    echo "       The sampling cache key omits the mesh, so two rows sharing an L"
    echo "       would reuse each other's samples. Give each row its own L."
    exit 1
fi

logdir="$outdir/logs"
mkdir -p "$outdir" "$logdir"

echo "== Exp4: target ${TARGET}^2, $NLEVELS levels, outdir=$outdir no_compute=$no_compute =="
for i in "${!MX_values[@]}"; do
  mx=${MX_values[$i]}; J=${J_values[$i]}; L=${L_values[$i]}
  # Define output file name based on the current h_0 value.
  output_file="mx_val_${mx}"
  # Define a clear log file
  log_file="$logdir/hLfixedNlevelsFixedMeshx${mx}"
  echo "-- mesh_x=mesh_y=$mx  J=$J  L=$L  grids $(( mx * (1 << J) ))..$(( mx * (1 << L) ))  -> $log_file"
  # Execute the python command with the current mesh_x value and redirect the output to the corresponding log file.
  python driver_MLCSPG_modulated_2D.py --cfg ../data/Exp4_meshxvaries.ini \
      --experiment_name "$outdir" --no_compute $no_compute \
      --mesh_x "$mx" --mesh_y "$mx" \
      --nb_level "$L" --l_start "$J" --nb_tests 0 \
      --output_file "$output_file" > "$log_file"
done
