#!/bin/bash
#
# launchGT_DRAFT.sh -- generate the reference ("ground truth") samples on a
# SLURM cluster.  Site-neutral template: fill the PLACEHOLDERS block below.
#
# Three stages, chained by afterok, so a failure stops the chain:
#   1. calibrate  : one short job. Warms the FFC/dijitso JIT cache AND runs
#                   compute_ground_truth.py --calibrate to MEASURE the solver
#                   tolerance, the preconditioner reuse, s/solve and the memory
#                   footprint, instead of guessing them.
#   2. shards     : an array, one task per shard of the sample set. The samples
#                   are independent, so this is the only real speed lever.
#   3. merge      : --merge, once every shard is done.
#
#   bash launchGT_DRAFT.sh --dry-run        # print the sbatch commands only
#   bash launchGT_DRAFT.sh --calibrate-only # stage 1 alone: do this FIRST
#   bash launchGT_DRAFT.sh                  # all three stages
#
# STRATEGY picks which of the two generators runs in stages 2 and 3:
#   single : compute_ground_truth.py on ONE mesh (MESH)
#   rich   : richardson_extrapolation.py on the MESHES ladder, which reaches
#            O(h^3) from meshes that a direct solver can still handle, and which
#            reports how far its finest mesh sits from the extrapolated limit
#            (i.e. it measures the quality a 'single' run would have to assume)
# Stage 1 is compute_ground_truth.py --calibrate either way, on the coarsest
# mesh in play: that is the only mesh where a direct solve is cheap enough to
# calibrate the iterative tolerance against, and it warms the JIT cache.
#
# Re-running is safe and cheap: a run reads <prefix>yGT.npy, works out what is
# missing, and shares only that out over the shards. A timed-out shard, a
# deleted shard file or a larger --nb_samples all just mean less or more work,
# never wrong values.
#
# NOTE: this script is young and has not been reviewed line by line. Run
# --dry-run and --calibrate-only before trusting it with a long queue.
# ---------------------------------------------------------------------------

set -u

# ============================= PLACEHOLDERS =================================
# Everything here is site-specific. The script REFUSES TO RUN while any
# <<<SET-ME>>> remains and tells you which ones are left.
#
#   PARTITION     SLURM partition(s), comma-separated.
#                   sinfo -s                       # names, nodes, TIMELIMIT
#                 Pick one whose TIMELIMIT covers TIME_SHARD below. Listing
#                 several ("short,long") lets SLURM take whichever frees first.
#
#   CONDA_BASE    the conda/mamba prefix, i.e. the directory containing
#                 etc/profile.d/conda.sh.
#                   conda info --base
#                 Needed because a batch job is non-interactive, so ~/.bashrc
#                 typically returns before conda reaches PATH. Can also be
#                 passed as an environment variable instead of edited here.
#
#   ACCOUNT       the accounting project to charge, "" if your site has none.
#                   sacctmgr show assoc user=$USER format=account
#
#   EXTRA_SBATCH  any further site flags, as a bash array. Examples:
#                   EXTRA_SBATCH=(--qos=normal)
#                   EXTRA_SBATCH=(--constraint=cascadelake --exclusive)
#                 Leave as () if you need none.
PARTITION="<<<SET-ME>>>"
CONDA_BASE="${CONDA_BASE:-<<<SET-ME>>>}"
ACCOUNT=""
EXTRA_SBATCH=()
# ============================================================================

# ========================= problem settings =================================
# These describe WHAT to compute. They are not site-specific, but every one of
# them is part of the reference's identity (see utilities/gt_meta.py), so
# changing any of them produces a different output folder -- by design.
STRATEGY=rich          # single | rich   (see the header)
MESH=3000              # 'single': the reference mesh width (both directions)
MESHES="640 1280 2560" # 'rich':   the ladder, coarsest first, each twice the previous
NB_SAMPLES=2000        # enough for the MEAN error; keep it fixed across curves
SEED=20261106          # MUST stay fixed: every curve must use the same points
CFG="../data/Exp1_hfvaries.ini"   # the coefficient, QoI and degree come from HERE
SOLVER=auto            # auto | cholmod | cg-amg | umfpack | ...
GT="compute_ground_truth.py"
RICH="richardson_extrapolation.py"
ENVNAME="env-fenics38" # the conda env; created by ../../prep_env_fenics38.sh
# ============================================================================

# ===================== resources: MEASURE, do not guess =====================
# Run with --calibrate-only first. Stage 1 writes a report to
# calibration/launchGT_<cfg>_h<mesh>/ comparing cg-amg against a DIRECT solve
# on the same sample points, and the numbers below come straight off it.
#
#   RTOL        from table 1. Pick the loosest tolerance whose max |dG| sits
#               at least two orders of magnitude below the smallest error you
#               intend to report. Accuracy is usually nearly free here: the AMG
#               setup dominates a solve, not the CG loop, so a much tighter
#               rtol often costs only a few percent.
#   PC_REBUILD  from table 2. The largest k whose max |dG| still matches the
#               k=1 row and that never prints DIVERGED. 1 = rebuild the AMG
#               hierarchy every sample (always safe), 0 = build it once.
#               Expect a modest win (tens of percent), not a factor.
#   TIME_SHARD  s/solve from the report, x NB_SAMPLES/NSHARDS, scaled from the
#               calibration mesh to the real one by the ratio of dof counts
#               (cost ~ dofs for an iterative solver), then x2-4 for headroom.
#   MEM_SHARD   the report prints the affine cache size; scale it the same way
#               and add the AMG hierarchy and the matrix values. Verify after
#               the first shard with:
#                 sacct -X -o JobName%10,Elapsed,MaxRSS --units=G -S today
#   NSHARDS     free choice: more shards = shorter wall time, same CPU total.
#               Keep NB_SAMPLES/NSHARDS well above 1.
#   CPUS        1 is right for a shard: it only solves PDEs, and nothing in
#               that path threads usefully. (Sparse recovery does, but that is
#               the experiment driver's job, not this script's.)
#
# The values below are deliberately generous placeholders, NOT measurements.
RTOL=1e-12
PC_REBUILD=1
NSHARDS=16
CPUS=1
MEM_SHARD="32G"
TIME_SHARD="06:00:00"
TIME_CAL="06:00:00"    # one direct row + 11 iterative rows of CAL_SAMPLES solves
CAL_MESH=1280          # stage 1 runs HERE, not on MESH: it needs a DIRECT solve
                       # to compare against, which stops fitting on a fine mesh.
                       # Lower it if stage 1 is slow; raise it to calibrate
                       # closer to the real mesh. It is capped by the coarsest
                       # mesh actually in play.
CAL_SAMPLES=20         # solves per calibration row; 10 is usually plenty
# ============================================================================

check_placeholders() {
    __bad=0
    for __v in PARTITION CONDA_BASE; do
        if [[ "${!__v}" == *"<<<SET-ME>>>"* ]]; then
            echo "FATAL: $__v is still <<<SET-ME>>>. See the PLACEHOLDERS block"
            echo "       at the top of $(basename "$0") for how to find its value."
            __bad=1
        fi
    done
    (( __bad == 0 )) || exit 2
}

# Site flags that every stage shares. Built as an array so empty values vanish
# instead of turning into an empty argument.
sbatch_site() {
    printf '%s\n' --partition="$PARTITION" ${ACCOUNT:+--account="$ACCOUNT"} \
                  ${EXTRA_SBATCH[@]+"${EXTRA_SBATCH[@]}"}
}

# ---- solver choice: by mesh only, deliberately -----------------------------
# A direct factorisation is exact but its memory grows faster than the problem:
# compute_ground_truth.report_size models a CHOLMOD factor at ~2.25 kB/dof with
# a ~1.9x peak during factorisation (mesh 1280 -> ~7 GB, 2560 -> ~32 GB). Above
# the crossover the iterative solver is the only option, and --calibrate is
# what tells you its error.  Set SOLVER explicitly to override.
pick_solver() {
    if [[ "$SOLVER" != "auto" ]]; then echo "$SOLVER"; return; fi
    if (( $(finest_mesh) <= 1280 )); then echo "cholmod"; else echo "cg-amg"; fi
}

# The two meshes the rest of the script needs: the finest decides the solver and
# the memory, the coarsest is where stage 1 calibrates.
finest_mesh() {
    if [[ "$STRATEGY" == "rich" ]]; then echo "${MESHES##* }"; else echo "$MESH"; fi
}
coarsest_mesh() {
    if [[ "$STRATEGY" == "rich" ]]; then echo "${MESHES%% *}"; else echo "$MESH"; fi
}
cal_mesh() {
    __c=$(coarsest_mesh)
    if (( __c < CAL_MESH )); then echo "$__c"; else echo "$CAL_MESH"; fi
}
# Stage 1's output goes here, deliberately nowhere near the GTresults_* folders:
# these are advisory numbers used to pick RTOL and PC_REBUILD, not reference
# values. The report inside says so in a banner, and carries the identity hash
# in its file name so two coefficient sets cannot overwrite each other.
cal_dir() {
    __cfg=$(basename "$CFG"); __cfg=${__cfg%.ini}
    echo "calibration/launchGT_${__cfg}_h$(cal_mesh)"
}
gt_script() {
    if [[ "$STRATEGY" == "rich" ]]; then echo "$RICH"; else echo "$GT"; fi
}
mesh_args() {
    if [[ "$STRATEGY" == "rich" ]]; then echo --meshes $MESHES; else echo --mesh "$MESH"; fi
}

# mpirun binds a single rank to ONE core by default (OpenMPI --bind-to core),
# which pins every numba and BLAS thread to that same core. Harmless for a
# single-threaded shard, fatal for anything threaded, so turn it off whenever
# this MPI accepts the flag. Delete this and the $MPIOPTS uses below if your
# site launches jobs with srun instead.
mpi_opts() {
    for opt in "--bind-to none" "-bind-to none"; do
        if mpirun $opt -n 1 true >/dev/null 2>&1; then echo "$opt"; return; fi
    done
}

activate_env() {
    __base=$(conda info --base 2>/dev/null)
    [[ -r "${__base:-}/etc/profile.d/conda.sh" ]] || __base="$CONDA_BASE"
    [[ -r "${__base}/etc/profile.d/conda.sh" ]] \
        || { echo "FATAL: no conda.sh under '${__base}'. Set CONDA_BASE to the"
             echo "       output of 'conda info --base'."; exit 1; }
    echo "### conda base: ${__base}"
    # conda's own activate.d/deactivate.d hooks read unset variables (e.g.
    # CONDA_BACKUP_CXX), which is fatal under 'set -u'. Relax it across the
    # activation only.
    set +u
    source "${__base}/etc/profile.d/conda.sh"
    conda activate "$ENVNAME"; __rc=$?
    set -u
    (( __rc == 0 )) || { echo "FATAL: conda activate $ENVNAME failed"; exit 1; }
    echo "### python: $(command -v python)"
}

common_args() {
    # no --out: the script derives GTresults_d<d>H<mesh>_<hash> itself, so a
    # different coefficient can never land in this one's folder
    echo --cfg "$CFG" $(mesh_args) --nb_samples "$NB_SAMPLES" \
         --seed "$SEED" --solver "$(pick_solver)" \
         --rtol "$RTOL" --pc_rebuild_every "$PC_REBUILD"
}

# ================================================ job modes (via sbatch) ===
if [[ "${1:-}" == "--stage" ]]; then
    stage="$2"
    check_placeholders
    # sbatch runs a COPY of this file from the spool directory, so $0 is useless
    # here and every relative path (CFG, the drivers, the output) needs this cd.
    cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
    activate_env
    # A private JIT cache: concurrent tasks otherwise race inside dijitso
    # ("RuntimeError: Destination file should exist at this point!").
    export DIJITSO_CACHE_DIR="$PWD/.dijitso/gt"
    mkdir -p "$DIJITSO_CACHE_DIR"
    MPIOPTS=$(mpi_opts)
    echo "### stage=$stage strategy=$STRATEGY mesh=$(mesh_args) solver=$(pick_solver) cpus=${SLURM_CPUS_PER_TASK:-?} host=$(hostname)"
    echo "### usable cores: $(mpirun $MPIOPTS -n 1 python -c \
        'import os; print(len(os.sched_getaffinity(0)))' 2>/dev/null)"

    case "$stage" in
      calibrate) mpirun $MPIOPTS -n 1 python "$GT" --cfg "$CFG" \
                     --mesh "$(cal_mesh)" --calibrate_samples "$CAL_SAMPLES" \
                     --seed "$SEED" --rtol "$RTOL" \
                     --calibrate_out "$(cal_dir)" --calibrate ;;
      shard)     mpirun $MPIOPTS -n 1 python "$(gt_script)" $(common_args) \
                     --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" ;;
      merge)     mpirun $MPIOPTS -n 1 python "$(gt_script)" $(common_args) \
                     --nshards "$NSHARDS" --merge ;;
      *) echo "unknown stage '$stage'"; exit 2 ;;
    esac
    exit $?
fi

# ====================================================== launcher mode ======
DRY=0; CALONLY=0
for a in "$@"; do
    case "$a" in
        --dry-run)        DRY=1 ;;
        --calibrate-only) CALONLY=1 ;;
        *) echo "unknown option: $a"; exit 2 ;;
    esac
done
(( DRY )) || check_placeholders

mkdir -p logs "$(cal_dir)"
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
run() { if (( DRY )); then printf '  %s' "$@"; echo; else "$@"; fi; }

case "$STRATEGY" in single|rich) ;; *) echo "STRATEGY must be single or rich"; exit 2 ;; esac
echo "== $STRATEGY: $(gt_script) $(mesh_args) samples=$NB_SAMPLES shards=$NSHARDS solver=$(pick_solver) cfg=$CFG =="
echo "   (the output folder is derived from the identity hash; each log prints it)"

echo "== 1/3  calibrate on mesh $(cal_mesh) -> $(cal_dir)/ (tolerance + preconditioner reuse) =="
mapfile -t SITE < <(sbatch_site)
cal=(sbatch --parsable --job-name=GTcal "${SITE[@]}"
     --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --time="$TIME_CAL" --mem="$MEM_SHARD"
     --output="$(cal_dir)/GTcal_%j.out" --error="$(cal_dir)/GTcal_%j.err"
     "$SELF" --stage calibrate)
if (( DRY )); then run "${cal[@]}"; CAL="<calibrate-jobid>"
else CAL=$("${cal[@]}") || exit 1; echo "  calibrate job id = $CAL"; fi

if (( CALONLY )); then
    echo
    echo "Read $(cal_dir)/GTcal_${CAL}.out (the same tables are kept as"
    echo "$(cal_dir)/calibration_<hash>.txt), then set RTOL, PC_REBUILD,"
    echo "TIME_SHARD and MEM_SHARD at the top of this script and run it again"
    echo "without --calibrate-only."
    exit 0
fi

echo "== 2/3  $NSHARDS shards, afterok:$CAL =="
shd=(sbatch --parsable --job-name=GTshard "${SITE[@]}"
     --nodes=1 --ntasks=1 --cpus-per-task="$CPUS"
     --time="$TIME_SHARD" --mem="$MEM_SHARD"
     --array="0-$((NSHARDS - 1))" --dependency="afterok:$CAL"
     --output=logs/GTshard_%A_%a.out --error=logs/GTshard_%A_%a.err
     "$SELF" --stage shard)
if (( DRY )); then run "${shd[@]}"; ARR="<shard-array-jobid>"
else ARR=$("${shd[@]}") || exit 1; echo "  shard array id = $ARR"; fi

echo "== 3/3  merge, afterok:$ARR =="
run sbatch --parsable --job-name=GTmerge "${SITE[@]}" \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:20:00 --mem=16G \
    --dependency="afterok:$ARR" \
    --output=logs/GTmerge_%j.out --error=logs/GTmerge_%j.err \
    "$SELF" --stage merge

echo
echo "watch:    squeue -u \$USER"
echo "progress: grep -h 'eta' logs/GTshard_*_*.out | tail -20"
echo "usage:    sacct -X -o JobID,JobName%10,Elapsed,MaxRSS,State --units=G -S today"
