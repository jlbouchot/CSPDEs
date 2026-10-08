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
    └── utilities/                          # parameter parsing, timing, ground truth, graphs
```


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
conda activate env-fenics38
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
a finer mesh, cached in `GTresults_d<2*nb_cosines>H<grid>/`.


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
