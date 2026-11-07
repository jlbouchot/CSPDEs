CSPDEs
======

Compressed sensing based approximation of solutions to high-dimensional parametric operator equations,
with a focus on the **Multi-Level Compressed Sensing Petrov-Galerkin (MLCSPG)** method.

> This toolbox is research code, provided as is. There is absolutely no guarantee that it will work on your problem!
>
> It relies on (legacy) FEniCS as a black box PDE solver. We cannot provide assistance with installing FEniCS,
> but we will happily help with questions about this implementation, ellipticity conditions, or compressed sensing.


About this project
==================

This project stems from the collaboration between
* Dr. Jean-Luc Bouchot, HDR, INRIA Saclay, GammaO Team (jean-luc.bouchot@inria.fr; previously at RWTH Aachen Chair C for Mathematics, then at Beijing Institute of Technology School of Mathematics and Statistics)
* Benjamin Bykowski, RWTH Aachen
* Prof. Holger Rauhut, Ludwig Maximilian University (then at Chair C for Mathematics (Analysis), RWTH Aachen)
* Prof. Christoph Schwab, Seminar for Applied Mathematics, ETH Zurich

What it does, in short: for a parametric diffusion problem whose coefficient depends on `d` parameters
(an affine cosine expansion with polynomially decaying weights), the code approximates a quantity of interest
(by default the average of the solution) as a sparse expansion in tensorized Chebyshev polynomials.
The coefficients are recovered from a few PDE solves at random parameter values, using weighted sparse recovery
algorithms (weighted IHT/HTP/OMP, or weighted basis pursuit through cvxpy). In the multi-level version, coarse
meshes are used with many samples and finer meshes with fewer samples, only to recover level-wise details.


Repository layout
=================

**Soon to be deprecated** Original version of the repo layout -- changes are in progress, see next paragraph

```
CSPDEs/
├── README.md, TODO.md
├── env-fenics38.yml                        # conda environment, Python 3.8 (recommended, see Installation)
├── env-fenics38.lock.txt                   # exact verified environment (linux-64)
├── prep_env_fenics38.sh                    # creates and tests env-fenics38
├── check_env_fenics38.py                   # environment check (imports, FEniCS solvers, cvxpy)
├── env-fenics.yml, working_conda_env.yml   # older Python 3.6 environments
├── prep_install_miniconda.sh               # Miniforge install (e.g. on a cluster, in $SCRATCH)
├── prep_env_fenics.sh                      # older Python 3.6 environment
├── notes_fenics-user                       # detailed installation notes and pitfalls
├── Papers/                                 # accompanying documents
└── CSPG/
    ├── CSPDE.py, Check.py                  # single-level algorithm and test harness
    ├── CSPDE_ML.py, Check_ML.py            # multi-level algorithm and test harness (writes results)
    ├── SPDE/                               # PDE side: FEniCS models, coefficients, goal functionals
    │   └── FEniCSModels/
    ├── WR/                                 # weighted recovery side
    │   ├── Algorithms/                     # wiht, whtp, womp, weighted BP / BPDN (cvxpy), ...
    │   ├── Operators/                      # Chebyshev (dense and tensor-based), Legendre, ...
    │   └── Weights/
    ├── data/                               # configuration files (*.ini)
    │   ├── mlcspg-default-cfg.ini          # defaults, always read first
    │   └── Exp*_*.ini                      # one file per experiment of the MLCSPG paper
    ├── testing/                            # drivers and launch scripts
    │   ├── driver_MLCSPG_2D.py             # main entry point (2D diffusion problem)
    │   ├── Exp*NoBatch.sh                  # run an experiment locally
    │   └── Exp*Batch.sh                    # same, as a SLURM array job
    ├── tests/                              # unit and regression tests (see tests/README.md)
    └── utilities/                          # parameter parsing, timing, ground truth, graphs
```

**New version** of the computations, with drastically increased speeds 
```
CSPDEs/
├── README.md, TODO.md
├── env-fenics38.yml                        # conda environment, Python 3.8 (recommended, see Installation)
├── env-fenics38.lock.txt                   # exact verified environment (linux-64)
├── prep_env_fenics38.sh                    # creates and tests env-fenics38
├── check_env_fenics38.py                   # environment check (imports, FEniCS solvers, cvxpy)
├── env-fenics.yml, working_conda_env.yml   # older Python 3.6 environments
├── prep_install_miniconda.sh               # Miniforge install (e.g. on a cluster, in $SCRATCH)
├── prep_env_fenics.sh                      # older Python 3.6 environment
├── notes_fenics-user                       # detailed installation notes and pitfalls
├── Papers/                                 # accompanying documents
└── CSPG/
    ├── CSPDE.py, Check.py                  # single-level algorithm and test harness
    ├── CSPDE_ML.py, Check_ML.py            # multi-level algorithm and test harness (writes results)
    ├── SPDE/                               # PDE side -- MAINLY DEPRECATED, changes in progress
    │   └── FEniCSModels/
    │       ├── *FEMModel*.py               # mesh + weak form + solver + QoI in ONE class
    │       └── WeightedCosSin2D.py, ...    # coefficient expansions -- still used below
    ├── SPDESampler/                        # replaces SPDE/: the same job split in four, which
    │                                       # is what buys the 20-40x (see the note below)
    │   ├── Discretisation/                 # mesh, function space, boundary conditions
    │   ├── Functionals/                    # goal functionals: AverageQOI, BoxAverage (--qoi)
    │   ├── Physics/                        # the PDE in weak form: Diffusion, NonlinearDiffusion
    │   ├── Samplers/                       # u(z) for one z; caches whatever does not depend on z
    │   └── Solvers/                        # one interface over umfpack / mumps / cholmod /
    │                                       # superlu / petsc-lu / cg-amg (iterative) / scipy
    ├── WR/                                 # weighted recovery side
    │   ├── Algorithms/                     # wiht, whtp, womp, weighted BP / BPDN (cvxpy), ...
    │   ├── Operators/                      # Chebyshev (dense and tensor-based), Legendre, ...
    │   └── Weights/
    ├── data/                               # configuration files (*.ini)
    │   ├── mlcspg-default-cfg.ini          # defaults, always read first
    │   └── Exp*_*.ini                      # one file per experiment of the MLCSPG paper
    ├── testing/                            # drivers and launch scripts
    │   ├── driver_MLCSPG_2D.py             # main entry point (2D diffusion problem)
    │   ├── Exp*NoBatch.sh                  # run an experiment locally
    │   └── Exp*Batch.sh                    # same, as a SLURM array job
    ├── tests/                              # unit and regression tests (see tests/README.md)
    └── utilities/                          # parameter parsing, timing, ground truth, graphs
```

Why the split in `SPDESampler/` is faster: nothing but the coefficient values depends on the parameter `z`.
The mesh, the function space, the sparsity pattern and the goal functional do not. Separating them lets each
be built **once** instead of once per sample:

* `Discretisation` is constructed once, so the mesh and the dof map are never rebuilt;
* `Solvers` keep their symbolic factorisation across samples (the pattern is fixed), so a sample pays only a
  numeric factorisation -- or, with `cg-amg`, only an AMG setup, which `--pc_rebuild_every` can reuse too;
* `AffineLinearSampler` caches the `d+1` matrices of the affine expansion, so assembling `A(z)` is a weighted
  sum of arrays rather than a FEniCS assembly. That cache is the memory driver at a fine mesh
  (`8*(d+1)*nnz` bytes, hence `--max_cache_gb`); `AssembledLinearSampler` reassembles instead when it does
  not fit.

The old `SPDE/FEniCSModels/` classes did all four jobs in one object and therefore redid all of it per
sample. Only the coefficient expansions (`WeightedCosSin2D`, `WeightedCosine2D`, ...) are shared between the
two layers, which is why `SPDE/` is still here.


Installation
============

Requirements
------------
* Linux. Tested on Ubuntu, Fedora 43 and on INRIA's HPC cluster [Margaret](https://clusters-saclay.gitlabpages.inria.fr/clusters-docs/docs/margaret/home/).
  Some tests have been done on LMU's supercomputing centre; reach out if you need to know more.
* [FEniCS](https://fenicsproject.org/) **2019.1.0 (legacy dolfin)**, used as the black box PDE solver.
  FEniCSx / dolfinx is **not** supported yet: this project started before it was released.
* [cvxpy](https://www.cvxpy.org/) 1.1, only for the convex recovery algorithms (`bp`, `bpdn`).
* numpy, progressbar2, matplotlib (graphs), pandas (some test scripts), numba (optional).

Since FEniCS 2019 is no longer maintained, recent versions of some packages are incompatible with it.
Pin versions as below and use strict channel priority.

Recommended: conda environment with Python 3.8
----------------------------------------------
If you do not have conda yet, `prep_install_miniconda.sh` installs Miniforge (it assumes `$SCRATCH` is set, as on a cluster).

```bash
bash prep_env_fenics38.sh            # creates env-fenics38 from env-fenics38.yml and tests it
conda activate env-fenics38
```

Or by hand:
```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f env-fenics38.yml
python check_env_fenics38.py         # imports, FEniCS solves with several linear solvers, cvxpy
```

`env-fenics38.yml` uses conda-forge only, with versions pinned to a verified environment.
Do not mix in packages from the `defaults` channel. To rebuild the exact verified environment
(linux-64 only), run `bash prep_env_fenics38.sh --exact`, which uses `env-fenics38.lock.txt`.

Try the following sequence of commands: 

Getting the code
----------------
Nothing to build, simply clone the repository:
```bash
git clone https://github.com/jlbouchot/CSPDEs.git
cd CSPDEs/CSPG
```


Quick start
===========

All drivers and launch scripts are meant to be run **from `CSPG/testing/`** (paths to configs and outputs are relative).

Smallest end-to-end check, using the debug sizes of experiment 3:
```bash
cd CSPG/testing
conda activate env-fenics38                # Or whatever environment works for you
DEBUG_MODE=true bash Exp3NoBatch.sh        # results and logs go to results/Exp3_debug/
```

Only check the sizes involved (sparsity per level, number of samples, size of the Ansatz space), without any PDE solve:
```bash
DEBUG_MODE=true NO_COMPUTE=true bash Exp3NoBatch.sh   # -> results/Exp3_debug_nocompute/
```

Run the driver directly:
```bash
python driver_MLCSPG_2D.py --cfg ../data/Exp3_Jvaries.ini \
    --experiment_name results/my_test --output_file my_run \
    --nb_level 4 --l_start 2 -e 0.000125
```

`bash launchAll.sh` runs the debug version of experiments 1 to 3 in sequence.

When running the conda environment non-interactively, prefer `conda run --no-capture-output -n env-fenics38 ...`
so that the progress output is not buffered.


Configuration
=============

Parameters are collected in three layers, each overriding the previous one:
1. `CSPG/data/mlcspg-default-cfg.ini` (always read),
2. the experiment file given with `--cfg`,
3. command line options (`python driver_MLCSPG_2D.py --help` lists them all).

The `.ini` files have three sections:

| Section          | Main keys | Meaning |
|------------------|-----------|---------|
| `[main]`         | `nb_cosines` (d), `power`, `abar`, `fluctuation_importance`, `weight_cosine` | parametric diffusion coefficient |
|                  | `n`, `mesh_x`, `mesh_y` | dimension and number of grid points of the coarsest mesh |
|                  | `nb_level` (L), `l_start` (J) | finest level, and first level used |
|                  | `gamma`, `exponent` | weights v_j = gamma * j^exponent |
|                  | `dat_constant`, `const_sj`, `p_0`, `p_t`, `t_0`, `t_prime` | sparsity per level s_l and the smoothness assumptions it is derived from |
|                  | `sampling` (`new`, `theoretic`, `pragmatic`) | rule for the number of samples m from s and the Ansatz space size N |
|                  | `do_tensor` | tensor-based (matrix free) Chebyshev operator |
|                  | `experiment_name`, `output_file` | output folder and result name (see below) |
|                  | `no_compute` | only print the sizes, no PDE solve nor recovery |
| `[pdesolver]`    | `linear_solver`, `preconditioner`, `elements`, `degree` | FEniCS discretization and linear solver |
| `[sparsesolver]` | `recovery_algo` (`wiht`, `whtp`, `womp`, `bp`, `bpdn`), `nb_iter`, `tol_res` | weighted sparse recovery |

Note that `tol_res` (CLI: `-e`) is used as the stopping criterion of the iterative recovery algorithms.
If it is larger than the size of the details on a level, the algorithm stops immediately and that level contributes nothing.

Every run writes the configuration it actually used to `<output_file>_config_file.txt`, next to its results.


Outputs
=======

Everything for one experiment goes to the folder `experiment_name` (relative to where the driver is run):

| File | Content |
|------|---------|
| `<output_file>` (`.dat/.dir/.bak` or `.db`, depending on the system) | Python `shelve` with the `TestResult` (models, L, per-level recovery results) |
| `<output_file>_config_file.txt` | the three configuration dictionaries, as JSON |
| `sampling_points_*.npy`, `sampling_points_*_y.npz` | cached sample points and PDE solves, per level and sparsity |
| `*_univariate_datamtx_tensor_*.npy` | cached tensor-based measurement matrices |
| `logs/` | stdout of each run (when using the `Exp*` scripts) |

The cached samples are **reused** by later runs in the same folder with the same `d`, level and sparsity,
which avoids recomputing expensive PDE solves. Use a new `experiment_name` to start from scratch.

If results named `<output_file>` already exist, a warning is printed: the new run adds an entry to the existing
shelve, and the graphing scripts read all entries.

To read results back:
```python
import shelve
results = sorted(shelve.open("results/Exp3_debug/J_val_2").values(), key=lambda r: r.L)
```

Results, caches, logs and archives are excluded from git (see `.gitignore`).


Experiments (MLCSPG paper)
==========================

WARNING: This section is being reviewed. More updates coming soon.

| Experiment | Config | Varies | Local / SLURM scripts |
|------------|--------|--------|-----------------------|
| 1 | `Exp1_hfvaries.ini` | target accuracy (finest level), h_0 and number of levels fixed | `Exp1NoBatch.sh` / `Exp1Batch.sh` |
| 2 | `Exp2_h0varies.ini` | starting mesh h_0, target accuracy and number of levels fixed | `Exp2NoBatch.sh` / `Exp2Batch.sh` |
| 3 | `Exp3_Jvaries.ini`  | first level J, h_0 and target accuracy fixed | `Exp3NoBatch.sh` / `Exp3Batch.sh` |
| 4 | `Exp4_Jvaries.ini`  | first level J, finer coarsest mesh | `Exp4Batch.sh` |
| 5 | `Exp5_dvaries.ini`  | number of parameters d | `Exp5NoBatch.sh` / `Exp5Batch.sh` |

Each script accepts the environment variables:
* `DEBUG_MODE=true`: small sizes, outputs to `results/ExpX_debug`;
* `NO_COMPUTE=true`: only check the sizes, outputs to `results/ExpX[_debug]_nocompute`.

On a SLURM cluster, submit the batch version, e.g. `sbatch Exp3Batch.sh`: each value of the varying parameter
is one task of the array job. Make sure `#SBATCH --array` matches the number of values in the script.

Graphs are produced by the scripts in `CSPG/utilities/` (`graphWithRespectToJ.py`, `graphWithRespectToH0.py`,
`graphTargetL.py`, `graphMLvsSL.py`, `graphCheckingDim.py`). They compare against a ground truth computed on
a finer mesh — see the next section.


Ground truth
============

Every error curve is measured against reference values `y = QoI(u(z))` computed on a mesh far finer than any
level of the experiment. Two scripts produce them, both run from `CSPG/testing/`:

| Script | Strategy |
|--------|----------|
| `compute_ground_truth.py` | one fine mesh (default 3000x3000) |
| `richardson_extrapolation.py` | extrapolate in h from a ladder of coarser meshes (e.g. 640/1280/2560), which reaches O(h^3) with meshes a direct solver can still handle |

Both take the coefficient, the QoI and the element degree from the experiment's `.ini`:

```bash
cd CSPG/testing
python compute_ground_truth.py --cfg ../data/Exp1_hfvaries.ini --mesh 3000 --nb_samples 2000
bash launchGT.sh --calibrate-only     # on SLURM: measure s/solve and memory first
bash launchGT.sh                      # then calibrate -> sharded array -> merge
```

### Calibrating before the big run

Two settings cannot be chosen from theory, so stage 1 of `launchGT.sh` measures them against a **direct**
solve on the same sample points, on a mesh where that is still affordable (`CAL_MESH`, default 1280 — never
on `MESH` itself when that is 3000²):

1. `RTOL` — the QoI error is bounded by `||w|| ||r||`, not by `rtol` alone. With condition number growing
   like `h^-2`, `rtol=1e-12` at 3000² admits a worst-case QoI error around `1e-5`, which is *larger* than the
   errors being plotted. The table gives the actual `max |dG|` per tolerance.
2. `PC_REBUILD` (`--pc_rebuild_every`) — `cg-amg` rebuilds the BoomerAMG hierarchy on every sample, although
   the sparsity pattern never changes between samples, only the coefficient values. Reusing it is the main
   speed lever at a fine mesh; the cost is CG iterations, and past some point the tolerance itself. `1` =
   rebuild every sample (default), `k` = every k-th, `0` = build once.

```bash
cd CSPG/testing
bash launchGT.sh --calibrate-only        # one job, CAL_SAMPLES solves x 12 configurations
# read calibration/launchGT_<cfg>_h<mesh>/, set RTOL and PC_REBUILD in launchGT.sh
bash launchGT.sh                         # then the full chain
```

Stage 1 writes **only** to its own folder, well away from any `GTresults_*`:

```
CSPG/testing/calibration/launchGT_Exp1_hfvaries_h1280/
├── calibration_<hash>.txt    # the two tables, under a banner saying what this is
└── GTcal_<jobid>.out/.err    # the SLURM logs
```

The report opens with "CALIBRATION OUTPUT — there are NO reference values in this folder … safe to delete at
any time", followed by the identity it was calibrated *for*, so a folder left lying around needs no
remembering. The hash is in the file name only so two coefficient sets cannot overwrite each other. The whole
`calibration/` tree is gitignored.

Pick the loosest `rtol` whose `max |dG|` is two orders of magnitude below the smallest error you report, and
the largest `k` whose `max |dG|` still matches the `k=1` row at that tolerance and that never prints
`DIVERGED`. If no reuse row is both faster and as accurate, leave `PC_REBUILD=1`.


### A reference belongs to exactly one problem

`y` depends on far more than `d` and the mesh: `nb_cosines` (hence d), `power`, `weight_cosine`,
`fluctuation_importance`, `abar`, the right-hand side, the QoI, the mesh, the element degree and the sample
seed. Change any one of them and every stored value is meaningless — with no shape mismatch and no exception,
just wrong numbers. So the folder name carries an 8-character hash of all ten (`utilities/gt_meta.py`):

```
CSPG/testing/
└── GTresults_d20H3000_5cf055c9/          # d, mesh, hash of the full identity
    ├── GT_meta.json                      # the ten fields in full, + the convergence claim
    ├── GT_samples.npy                    # (nb_samples, d) the points, cos(pi U), seeded
    ├── GT_yGT.npy                        # (nb_samples,) the reference values  <- what the graphs read
    └── GT_yGT.shard<k>of<n>.npy          # scratch, one per shard; --merge folds them into GT_yGT.npy
```

`GT_yGT.npy` is always `nb_samples` long, with `NaN` where a sample is still missing, and it is the only
authoritative file: it is never shrunk, and the shard files are disposable scratch (each value is stored next
to the sample *index* it belongs to).

`--merge` **deletes the shard files** once it has folded them in (`--keep_shards` opts out). They are
redundant at that point whether or not the set is complete, and a shard that reruns afterwards derives its
work from `GT_yGT.npy`, so nothing is re-solved. Keeping them would also let a value computed at a looser
`--rtol` be re-injected after a deliberate recompute with a direct solver, since the solver is deliberately
*not* part of the identity.

`richardson_extrapolation.py` writes the same folder shape, with the ladder in place of the single width, plus
the per-mesh values so the extrapolation can be redone without re-solving:

```
└── GTresults_d20Hrich640-1280-2560_a1b2c3d4/
    ├── GT_meta.json                      # also records the MEASURED order (should be ~2)
    ├── GT_samples.npy
    ├── GT_levels.npz                     # (nb_meshes, nb_samples) the QoI on each mesh
    └── GT_yGT.npy                        # the extrapolated reference
```

A run into a folder whose `GT_meta.json` describes a different problem **aborts** and names the folder it
should have used.

### Extending a set

`nb_samples` is deliberately *not* part of the identity, and the draw is prefix-stable, so a set may simply
grow. Raise
`NB_SAMPLES` and rerun: each run reads `GT_yGT.npy`, works out what is missing (the new tail, plus any `NaN`
left by a failed solve), shares **only the missing samples** over the shards, and leaves everything else
alone. Going from 2000 to 3000 costs 1000 solves whether or not the shard files from the first round are
still there. Keep `MESH`, `SEED` and the `.ini` identical or the identity changes and you get a different
folder — which is the point.

A shard takes its slice from `GT_yGT.npy` *before* looking at any shard file, so all the tasks of one round
agree on the same partition even though SLURM starts them minutes apart. (An earlier version sliced the
missing set as each task saw it, which left samples that no shard owned.)


### Failures and restarts

A sample that does not converge is stored as `NaN` and the run **continues** — one bad solve must not throw
away the hours already spent on the others. At the end the script lists the offending indices and exits 1
(unless `--allow_failures`), so the `afterok` merge does not run on a set with holes.

There is **no automatic retry**: rerunning the same command is the retry. A restart recomputes exactly the
`NaN`s and nothing else, so after loosening `--rtol` or switching to a direct `--solver` a rerun fixes just
the bad samples. The same mechanism covers a shard that hits its SLURM time limit, and `--merge` is
idempotent — it can be rerun as more shards finish, and will never remove a value it already has.

### Reading it back

The graph scripts do not spell a folder name out. Each declares, near the top,

```python
GT_CFG  = .../data/Exp2_h0varies.ini      # which experiment, hence which reference
GT_SEED = 20261106
```

and `utilities/utils.loadGroundTruth` resolves the folder from that identity, verifies the hash, regenerates
the points from the seed and prints the stored solver/tolerance. It **only reads**: if the reference is
missing or too short it stops and prints the command to generate it. (It used to compute whatever was
missing on the spot with gmres at `rtol=1e-7`, three orders of magnitude looser than the errors being
plotted — silently, inside a plotting script.)


Accompanying papers - Theory
============================
The details of these methods can be found in the following publications:
* H. Rauhut and C. Schwab,
"Compressive sensing Petrov-Galerkin approximation of high-dimensional parametric operator equations",
Mathematics of Computation 86(304):661-700, 2017.
([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/csparampde.pdf))

* J.-L. Bouchot, B. Bykowski, H. Rauhut and C. Schwab,
"Compressed sensing Petrov-Galerkin approximations for parametric PDEs",
Sampling Theory and Applications 2015 (SampTA 15). ([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/SampTA15_BBRS.pdf))

* B. Bykowski,
"Weighted l1 minimization methods for numerical approximations of parametric PDEs under uncertainty quantification",
Master's Thesis, Chair C for Mathematics, RWTH Aachen, July 2015. ([Thesis](./Papers/bykowski_master.pdf))

* J.-L. Bouchot, H. Rauhut and C. Schwab,
"Multi-level Compressed Sensing Petrov-Galerkin discretization of high-dimensional parametric PDEs",
Submitted Jan. 2017. ([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/MLCSPG.pdf))

* J.-L. Bouchot,
"Weighted block compressed sensing for parametrized function approximation",
Submitted Nov. 2018. ([Preprint](https://arxiv.org/abs/1811.04598))


License
=======
More details very soon.


Tips
====
List the parameters available for the FEniCS solvers:
```python
from dolfin import *
info(LinearVariationalSolver.default_parameters(), True)
info(NonlinearVariationalSolver.default_parameters(), True)
```


Acknowledgments
===============
This work was partly supported by the European Research Council through the grant StG 258926. Part of this work was developed as J.-L. B. and H. R. were visiting the Hausdorff Research Center for Mathematics as part of the Hausdorff Trimester Program on Mathematics of Signal Processing. 
