#!/bin/bash
#
# launchExp1_DRAFT.sh -- submit the Exp1 (h_f varies) sweep on a SLURM cluster.
# Site-neutral template: fill the PLACEHOLDERS block below.
#
# Acts as BOTH the launcher (run it on the login node) and the job script
# (sbatch re-invokes this same file with --job / --warm).
#
#   bash launchExp1_DRAFT.sh --dry-run        # print the sbatch commands only
#   bash launchExp1_DRAFT.sh                  # submit the whole sweep
#   bash launchExp1_DRAFT.sh --debug          # 3 short jobs only
#   bash launchExp1_DRAFT.sh --only "6 7"     # resubmit a subset (see below)
#
# Why it is built this way
#   * one short WARM-UP job first: it JIT-compiles the UFL forms into a shared
#     FFC/dijitso cache. Without it, concurrent jobs race in
#     dijitso.system._lockfree_move_file -> "Destination file should exist at
#     this point!". Every real job waits on afterok of this one.
#   * each real job then COPIES the warm cache into a private directory, so the
#     jobs are fully independent (no race) and still never recompile.
#   * the warm-up uses (J,L)=(2,3), a two-level shape ON PURPOSE: its s/N/m
#     schedule hits a level with N<m, so it exercises the cs_ok=False
#     least-squares fallback that the real L-J=2 jobs never reach.
#   * one job per (J,L) instead of one array: the cost spans orders of
#     magnitude across L, and an array can only carry a single --time/--mem.
#   * --nb_tests 0: no cross-check here; the reference values are generated
#     separately and once, by launchGT_DRAFT.sh.
#   * tol_res is NOT passed: the code derives it per level from
#     utils.default_tolerance_sparse_recovery(L, L-J+1, mesh_x, t_0, t_prime).
#
# NOTE: this script is young and has not been reviewed line by line. Run
# --dry-run and --debug before trusting it with a long queue.
# ---------------------------------------------------------------------------

set -u

# ============================= PLACEHOLDERS =================================
# Everything here is site-specific. The script REFUSES TO RUN while any
# <<<SET-ME>>> remains and tells you which ones are left.
#
#   PARTITION     SLURM partition(s), comma-separated.
#                   sinfo -s                       # names, nodes, TIMELIMIT
#                 The largest TIME_FOR_L below must fit inside the partition's
#                 TIMELIMIT, or the job is rejected at submission.
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
#                   EXTRA_SBATCH=(--qos=long)
#                   EXTRA_SBATCH=(--constraint=cascadelake)
#                 Leave as () if you need none.
PARTITION="<<<SET-ME>>>"
CONDA_BASE="${CONDA_BASE:-<<<SET-ME>>>}"
ACCOUNT=""
EXTRA_SBATCH=()
# ============================================================================

# ========================= problem settings =================================
CFG="../data/Exp1_hfvaries.ini"
DRIVER="driver_MLCSPG_modulated_2D.py"   # the driver that honours numba_threads
OUTBASE="results/exp1_h0_nbDisc_fixed"   # one folder for the whole sweep; the
                                         # per-level caches live inside it
CACHE_ROOT="$PWD/.dijitso"
ENVNAME="env-fenics38"                   # created by ../../prep_env_fenics38.sh
# ============================================================================

# ================== threads: MEASURE on your own hardware ===================
# Only the sparse recovery threads usefully (through numba); the PDE sampling
# does not, whatever these are set to. So the question is just "where does
# numba stop scaling on this machine", and the answer is per-CPU.
#
# How to measure it, cheaply: run one short (J,L) with --numba_threads 1, 2, 4,
# 6 in four separate jobs with --cpus-per-task set to the largest, each writing
# to its OWN --experiment_name (a shared folder would make later runs reload
# the first one's cached samples AND its stored timings, so they would look
# identical). Then compare, per run, the wall time and the CPU/wall ratio on
# the "Recovery time" line. Take the knee: the largest count still near-linear.
#
# Two traps worth knowing:
#   * if every thread count gives the same wall time AND a CPU/wall ratio near
#     1, the process is pinned to one core -- see mpi_opts() below.
#   * raising BLAS threads can make things SLOWER through contention with numba.
#     Measure it rather than assuming more is better.
#
# The values below are a safe, conservative starting point, NOT a measurement.
CPUS=4                 # --cpus-per-task; must be >= max(NUMBA, BLAS, OMPT)
NUMBA=4                # numba_threads: the one that usually matters
BLAS=1                 # OPENBLAS/MKL threads
OMPT=1                 # OMP_NUM_THREADS
# ============================================================================

# ============== walltime and memory per L: MEASURE, do not guess =============
# The cost grows steeply and non-obviously with L: the finest level dominates,
# and in our runs sampling grew by roughly x6 per level. Projecting from small
# L underestimates badly -- we had L=6 and L=7 killed by limits derived that
# way, so the numbers below are deliberately generous.
#
# HOW TO SET THEM for your problem:
#   1. run --debug first (small L). Note Elapsed and MaxRSS:
#        sacct -X -o JobName%14,Elapsed,MaxRSS,State --units=G -S today
#   2. extrapolate the SAMPLING cost by the factor you observe between two
#      consecutive L (ours was ~6x), not by the total.
#   3. multiply by 2-4. A too-long limit only costs queue priority. A too-SHORT
#      one is expensive: the cache is per COMPLETED level, so a level killed in
#      the middle is redone from scratch -- and at the largest L that single
#      level is most of the job.
#   4. check the result fits your partition's TIMELIMIT (sinfo -s).
#
# IMPORTANT LIMITATION: there is no per-sample checkpoint inside a level. If the
# finest level of your largest L cannot finish inside the maximum walltime your
# site allows, no number of resubmissions will ever complete it; the sampling
# loop needs checkpointing first.
declare -A TIME_FOR_L MEM_FOR_L
#                        L = 2     3     4     5     6     7     8
TIME_FOR_L=( ["2"]="00:20:00" ["3"]="00:30:00" ["4"]="01:00:00" \
             ["5"]="02:00:00" ["6"]="12:00:00" ["7"]="2-00:00:00" ["8"]="7-00:00:00" )
MEM_FOR_L=(  ["2"]="8G"       ["3"]="8G"       ["4"]="8G"       \
             ["5"]="8G"       ["6"]="16G"      ["7"]="32G"      ["8"]="64G" )
# Resubmitting is safe: samples and matrices are cached per (Lmax, level, s, d)
# under $OUTBASE and a re-run resumes from the last COMPLETED level. The cache
# is keyed by Lmax, so nothing is shared between different L and each L pays for
# its own samples.
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

# Site flags that every job shares. Built as an array so empty values vanish
# instead of turning into an empty argument.
sbatch_site() {
    printf '%s\n' --partition="$PARTITION" ${ACCOUNT:+--account="$ACCOUNT"} \
                  ${EXTRA_SBATCH[@]+"${EXTRA_SBATCH[@]}"}
}

# mpirun binds a single rank to ONE core by default (OpenMPI --bind-to core),
# pinning every numba and BLAS thread to that same core: the thread settings
# above then have NO effect and a sweep over them looks completely flat. Find
# the spelling this MPI accepts and turn binding off. Delete this and the
# $MPIOPTS uses below if your site launches jobs with srun instead.
mpi_opts() {
    for opt in "--bind-to none" "-bind-to none"; do
        if mpirun $opt -n 1 true >/dev/null 2>&1; then echo "$opt"; return; fi
    done
}

# -------------------------------------------------- conda, usable everywhere
activate_env() {
    __base=$(conda info --base 2>/dev/null)
    if [[ -z "${__base:-}" || ! -r "${__base}/etc/profile.d/conda.sh" ]]; then
        source "$HOME/.bashrc" 2>/dev/null
        __base=$(conda info --base 2>/dev/null)
    fi
    if [[ -z "${__base:-}" || ! -r "${__base}/etc/profile.d/conda.sh" ]]; then
        __base="$CONDA_BASE"
    fi
    echo "### conda base: ${__base}"
    if [[ ! -r "${__base}/etc/profile.d/conda.sh" ]]; then
        echo "FATAL: no conda.sh under '${__base}'. Set CONDA_BASE at the top"
        echo "       of $(basename "$0") to the output of 'conda info --base'."
        exit 1
    fi
    # conda's own activate.d/deactivate.d hooks read unset variables (e.g.
    # CONDA_BACKUP_CXX in deactivate-gxx_linux-64.sh), which is fatal under
    # 'set -u'. Relax it across the activation only.
    set +u
    source "${__base}/etc/profile.d/conda.sh"
    conda activate "$ENVNAME"
    __rc=$?
    set -u
    (( __rc == 0 )) || { echo "FATAL: conda activate $ENVNAME failed ($__rc)"; exit 1; }
    python -c "import numpy, numba" \
        || { echo "FATAL: $ENVNAME activated but numpy/numba missing"; exit 1; }
    echo "### python: $(command -v python)"
}

# =================================== job modes (invoked by sbatch) =========
if [[ "${1:-}" == "--warm" ]]; then
    check_placeholders
    # sbatch runs a COPY of this file from the spool directory, so $0 is useless
    # here and every relative path (CFG, the driver, the output) needs this cd.
    cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
    CACHE_ROOT="$PWD/.dijitso"
    activate_env
    export DIJITSO_CACHE_DIR="$CACHE_ROOT/shared"
    mkdir -p "$DIJITSO_CACHE_DIR" "$OUTBASE"_jitwarm
    echo "### warm-up, DIJITSO_CACHE_DIR = $DIJITSO_CACHE_DIR"
    MPIOPTS=$(mpi_opts); echo "### mpirun options: '${MPIOPTS:-<default>}'"
    mpirun $MPIOPTS -n 1 python "$DRIVER" --cfg "$CFG" \
        --experiment_name "${OUTBASE}_jitwarm" --no_compute False \
        --nb_level 3 --l_start 2 -N 1 --nb_tests 0 \
        --numba_threads 1 --blas_threads 1 --output_file jitwarm
    # capture the status BEFORE the echoes, or the warm-up always reports
    # success and 7 jobs race against a cold cache
    rc=$?
    echo "### warm cache contents: $(find "$DIJITSO_CACHE_DIR" -type f | wc -l) files"
    echo "### warm-up exit code: $rc"
    exit $rc
fi

if [[ "${1:-}" == "--job" ]]; then
    J="$2"; L="$3"; no_compute="$4"; outdir="$5"
    check_placeholders
    cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
    CACHE_ROOT="$PWD/.dijitso"
    activate_env
    # private cache, seeded from the warm one: no race, no recompile
    export DIJITSO_CACHE_DIR="$CACHE_ROOT/L${L}"
    if [[ ! -d "$DIJITSO_CACHE_DIR" ]]; then
        cp -a "$CACHE_ROOT/shared" "$DIJITSO_CACHE_DIR" 2>/dev/null \
            || mkdir -p "$DIJITSO_CACHE_DIR"
    fi
    mkdir -p "$outdir"
    echo "### J=$J L=$L no_compute=$no_compute outdir=$outdir"
    echo "### cpus=${SLURM_CPUS_PER_TASK:-?} numba=$NUMBA blas=$BLAS omp=$OMPT host=$(hostname)"
    echo "### DIJITSO_CACHE_DIR = $DIJITSO_CACHE_DIR"
    MPIOPTS=$(mpi_opts); echo "### mpirun options: '${MPIOPTS:-<default>}'"
    # Sanity check worth reading in every log: if this prints 1 while CPUS > 1,
    # the thread settings above are doing nothing (see mpi_opts).
    echo "### usable cores: $(mpirun $MPIOPTS -n 1 python -c \
        'import os; print(len(os.sched_getaffinity(0)))' 2>/dev/null)"
    mpirun $MPIOPTS -n 1 python "$DRIVER" --cfg "$CFG" \
        --experiment_name "$outdir" --no_compute "$no_compute" \
        --nb_level "$L" --l_start "$J" --nb_tests 0 \
        --numba_threads "$NUMBA" --blas_threads "$BLAS" --omp_threads "$OMPT" \
        --output_file "L_val_${L}"
    exit $?
fi

# ======================================================= launcher mode =====
DRY=0
DEBUG=${DEBUG_MODE:+1}; DEBUG=${DEBUG:-0}
NOCOMP=${NO_COMPUTE:+1}; NOCOMP=${NOCOMP:-0}
ONLY=""
while (( $# )); do
    case "$1" in
        --dry-run)    DRY=1 ;;
        --debug)      DEBUG=1 ;;
        --no-compute) NOCOMP=1 ;;
        ## Resubmit a SUBSET, e.g. --only "6 7". Needed because every L shares
        ## $OUTBASE: two jobs for the same L would write the same
        ## sampling_points_Lmax<L>_* files at the same time. Different L values
        ## never collide, so running them side by side is fine.
        --only)       ONLY="$2"; shift ;;
        *) echo "unknown option: $1"; exit 2 ;;
    esac
    shift
done
(( DRY )) || check_placeholders

# outdir and no_compute are resolved HERE and handed to every job, so the
# launcher and the jobs can never disagree about where results land.
outdir="$OUTBASE"
no_compute=False
if (( DEBUG ));  then outdir="${outdir}_debug"; fi
if (( NOCOMP )); then outdir="${outdir}_nocompute"; no_compute=True; fi

# One (J,L) pair per job. L-J is held fixed at 2 here: the sweep is over the
# target accuracy (the finest level), with the number of levels fixed.
if (( DEBUG )); then
    J_values=(0 1 2); L_values=(2 3 4)
else
    J_values=(0 1 2 3 4 5 6); L_values=(2 3 4 5 6 7 8)
fi

mkdir -p logs "$CACHE_ROOT"
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
run() { if (( DRY )); then printf '  %s' "$@"; echo; else "$@"; fi; }
mapfile -t SITE < <(sbatch_site)

echo "== settings: outdir=$outdir no_compute=$no_compute cpus=$CPUS numba=$NUMBA blas=$BLAS omp=$OMPT =="
echo "== 1/2  warm-up job (fills $CACHE_ROOT/shared) =="
warm_cmd=(sbatch --parsable --job-name=Exp1warm "${SITE[@]}"
          --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:20:00 --mem=4G
          --output=logs/warm_%j.out --error=logs/warm_%j.err
          "$SELF" --warm)
if (( DRY )); then
    run "${warm_cmd[@]}"; WARM="<warmup-jobid>"
else
    WARM=$("${warm_cmd[@]}") || { echo "warm-up submit failed"; exit 1; }
    echo "  warm-up job id = $WARM"
fi

echo "== 2/2  one job per (J,L)${ONLY:+, restricted to L in [$ONLY]}, each afterok:$WARM =="
for i in "${!L_values[@]}"; do
    J=${J_values[$i]}; L=${L_values[$i]}
    if [[ -n "$ONLY" && " $ONLY " != *" $L "* ]]; then continue; fi
    t=${TIME_FOR_L[$L]:-01:00:00}; mem=${MEM_FOR_L[$L]:-8G}
    run sbatch --parsable --job-name="Exp1_L${L}" "${SITE[@]}" \
        --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" \
        --time="$t" --mem="$mem" --dependency="afterok:$WARM" \
        --output="logs/Exp1_L${L}_%j.out" --error="logs/Exp1_L${L}_%j.err" \
        "$SELF" --job "$J" "$L" "$no_compute" "$outdir"
done

echo
echo "watch:   squeue -u \$USER"
echo "timings: for f in logs/Exp1_L*_*.out; do echo \"== \$f\"; grep 'Recovery time' \$f; done"
echo "usage:   sacct -X -o JobID,JobName%14,Elapsed,MaxRSS,State --units=G -S today"
