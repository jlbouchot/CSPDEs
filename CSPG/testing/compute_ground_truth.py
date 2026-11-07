"""Precompute the reference ("ground truth") QoI samples.

Took the naming conventions from the utils.loadGroundTruth helper.
Writes exactly the two files that ``utilities.utils.loadGroundTruth`` looks for,

    <out>/<prefix>samples.npy     Z, shape (nb_samples, d), float64
    <out>/<prefix>yGT.npy         y, shape (k,), the first k reference values

so that the graph scripts can load the pre computed ground truth samples. 
The default sampling measure is reproduced exactly:
> Z = cos(pi U(-1,1)), i.e. LD_bounded_operator.apply_precondition_measure 
applied to a uniform draw, which is what the default usecases use in the accompanying paper. 

Usage
-----
Calibrate the tolerance first, on a mesh where the direct solver still fits:

    python compute_ground_truth.py --calibrate --mesh 1280 --nb_samples 8

then run the real thing, optionally as an array job:

    python compute_ground_truth.py --mesh 3000 --nb_samples 1000 --solver cg-amg --rtol 1e-12
    python compute_ground_truth.py --mesh 3000 --nb_samples 1000 --shard $i --nshards 16
    python compute_ground_truth.py --mesh 3000 --nb_samples 1000 --merge

Warning: This script was heavily vibe coded and not fully reviewed. It is a helper file, not a core part of the framework. It is not expected to be used outside of the experiments in the paper, and it is not expected to be maintained. Use at your own risk.
"""

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/11/06"
__lastmodified__ = "2026/11/07"

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np

from utilities import gt_meta


def add_common_options(p):
    """
    Common options shared between ground truth computation strategies. 
    Context: A ground truth can be computed at a high accuracy, e.g. 
    * via a direct solver
    * via a very tight tolerance iterative solver
    * via a Richardson extrapolation of two solves at different tolerances
    * ...
    The options here are options common to all gt computations. 
    Declared once so a new option cannot reach one script and not the other. This function
    needs to be called wherever a GT computation function is called. 
    """
    g = p.add_argument_group("reference discretisation")
    g.add_argument("--degree", type=int, default=1, help="Lagrange degree (default: 1)")

    g = p.add_argument_group("parametric coefficient (must match the experiment)")
    g.add_argument("--cfg", default=None,
                   help="experiment .ini to read the coefficient and QoI from "
                        "(nb_cosines, power, weight_cosine, fluctuation_importance, "
                        "abar, qoi, degree). Explicit command-line options win.")
    g.add_argument("--nb_cosines", type=int, default=16, help="number of frequencies; d = 2*nb_cosines")
    g.add_argument("--power", type=float, default=5.0, help="mu, the decay exponent")
    g.add_argument("--weight_cosine", type=float, default=0.5, help="Upsilon, the knee position")
    g.add_argument("--fluctuation_importance", type=float, default=3.9604513858, help="varsigma")
    g.add_argument("--abar", type=float, default=10.0, help="mean field")
    g.add_argument("--rhs", type=float, default=10.0, help="constant right-hand side f")
    g.add_argument("--qoi", choices=("box", "avg", "average"), default="box",
                   help="BoxAverage on [0,1/4]^2, or AverageQOI on the whole square")

    g = p.add_argument_group("samples")
    g.add_argument("--nb_samples", type=int, default=1000)
    g.add_argument("--seed", type=int, default=20261106,
                   help="fixes Z; every shard must use the same value")

    g = p.add_argument_group("solver")
    g.add_argument("--solver", default="cg-amg",
                   help="cg-amg (iterative), cholmod/mumps-cholesky/umfpack (direct)")
    g.add_argument("--rtol", type=float, default=1e-12, help="relative residual tolerance (iterative only)")
    g.add_argument("--atol", type=float, default=1e-16, help="absolute residual tolerance (iterative only)")
    g.add_argument("--assembly", choices=("auto", "affine", "per-sample"), default="auto",
                   help="'affine' caches the d+1 matrices and then costs one mat-vec per sample; "
                        "'per-sample' reassembles the form every sample, which at this mesh is "
                        "orders of magnitude slower. 'auto' takes affine whenever the cache fits "
                        "--max_cache_gb (default)")
    g.add_argument("--pc_rebuild_every", type=int, default=1,
                   help="cg-amg only: rebuild the AMG hierarchy every k-th sample "
                        "(1 = every sample, the default and PETSc's own behaviour; "
                        "0 = build it once). Reusing it is the main speed lever at a "
                        "fine mesh, but it costs CG iterations -- calibrate first")
    g.add_argument("--max_cache_gb", type=float, default=20.0,
                   help="budget for the affine cache, 8*(d+1)*nnz bytes (default: 20 GB of a 32 GB node)")

    g = p.add_argument_group("output and restart")
    g.add_argument("--out", default=None, help="default: GTresults_d<d>H<mesh>, as the graph scripts expect")
    g.add_argument("--prefix", default="GT_")
    g.add_argument("--allow_failures", action="store_true",
                   help="exit 0 even if some samples did not converge (default: "
                        "finish every sample, then exit 1 so an afterok merge "
                        "does not run on a set with holes)")
    g.add_argument("--checkpoint", type=int, default=25, help="flush to disk every this many samples")
    g.add_argument("--shard", type=int, default=None, help="this task's index, 0-based")
    g.add_argument("--nshards", type=int, default=None)
    g.add_argument("--merge", action="store_true", help="combine shard files into <prefix>yGT.npy and exit")
    g.add_argument("--keep_shards", action="store_true",
                   help="do not delete the shard files after --merge has folded "
                        "them in (they are redundant once it has: see merge())")
    return p


def parse_args_parser():
    """The parser alone, so gt_meta.check_parser can inspect it."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_options(p)
    # Part of the parser that is specific to a direct solver. 
    p.add_argument("--mesh", type=int, default=3000,
                   help="cells per direction of the reference mesh (default: 3000, as in the graph scripts)")
    p.add_argument("--calibrate", action="store_true",
                   help="compare the direct solver against cg-amg on --mesh (tolerances, then "
                        "preconditioner reuse), then exit")
    p.add_argument("--calibrate_out", default=None,
                   help="write the calibration report here as well as to stdout "
                        "(launchGT.sh passes calibration/launchGT_<cfg>_h<mesh>)")
    p.add_argument("--calibrate_samples", type=int, default=20,
                   help="samples per calibration row (default: 20). NOT --nb_samples: the "
                        "calibration solves this many times per configuration")
    return p


def parse_args(argv=None):
    return parse_args_parser().parse_args(argv)


def apply_cfg(args, argv):
    """Overlay args.cfg on the options the .ini also defines, then report where
    every value came from. Anything passed explicitly on the command line wins.
    The field list lives in gt_meta.CFG_FIELDS, next to the identity it feeds."""
    if not args.cfg:
        return args
    from utilities.MLCSPG_param_parser import (define_MLCSPG_cli_parser,
                                               parse_cfg_file_and_args)
    main_cfg, pde_cfg, _ = parse_cfg_file_and_args(
        define_MLCSPG_cli_parser().parse_args(["--cfg", args.cfg]))
    section = {"main": main_cfg, "pde": pde_cfg}
    given = set(a.split("=")[0] for a in (sys.argv[1:] if argv is None else argv))
    print("coefficient and QoI (cfg = %s):" % args.cfg)
    for k in sorted(gt_meta.CFG_FIELDS):
        cfg = section[gt_meta.CFG_FIELDS[k]]
        cur = getattr(args, k)
        if "--" + k in given:
            src = "command line (cfg says %s)" % cfg.get(k, "-")
        elif k in cfg:
            setattr(args, k, type(cur)(cfg[k]) if cur is not None else cfg[k])
            src = "cfg"
        else:
            src = "script default (absent from cfg)"
        print("  %-24s %-12s %s" % (k, getattr(args, k), src))
    print("  %-24s %-12s derived as 2 * nb_cosines" % ("d", 2 * args.nb_cosines))
    return args


def est_nnz(mesh, degree):
    """
    Free dofs and nonzeros of the P1 stiffness matrix, before building anything: the assembly
    strategy has to be chosen before the discretisation exists. Measured ratio 6.95 nnz/dof.
    """
    # TODO: Adapt in case of P2 or other stiffness matrices. 
    if degree != 1:
        return None, None
    free = (mesh - 1) ** 2
    return free, int(6.95 * free)


def choose_assembly(args, d):
    if args.assembly != "auto":
        return args.assembly
    free, nnz = est_nnz(args.mesh, args.degree)
    if nnz is None:
        return "per-sample"
    return "affine" if 8.0 * (d + 1) * nnz / 2**30 <= args.max_cache_gb else "per-sample"


def build(args, solver, assembly=None):
    """A sampler for the reference problem. ``solver`` is a name or a configured instance."""
    from SPDESampler.Discretisation import UnitDiscretisation
    from SPDESampler.Samplers import AffineLinearSampler, AssembledLinearSampler
    from utilities import utils as u
    from SPDESampler.Physics import Diffusion
    from SPDE.FEniCSModels.WeightedCosSin2D import WeightedCosSin2D

    d = 2 * args.nb_cosines
    coef = WeightedCosSin2D(d, args.power, args.fluctuation_importance, args.abar, args.weight_cosine)
    phys = Diffusion(coef, f=args.rhs)
    disc = UnitDiscretisation((args.mesh, args.mesh), degree=args.degree)
    qoi = u.get_functional_type(args.qoi)   # same mapping as the drivers
    assembly = assembly or choose_assembly(args, d)
    cls = AffineLinearSampler if assembly == "affine" else AssembledLinearSampler
    return cls(disc, phys, qoi, solver), disc, d


def make(args):
    """Instantiate the solver, passing the tolerances that BaseSampler's string form drops."""
    from SPDESampler.Solvers.linear_solvers import make_solver
    s = make_solver(args.solver, rtol=args.rtol, atol=args.atol,
                    pc_rebuild_every=args.pc_rebuild_every)
    return s


def report_size(disc, d, assembly):
    nnz = disc.indices.size
    free = disc.indptr.size - 1
    print("  reference mesh: %d free dofs, %d nonzeros" % (free, nnz))
    print("  matrix values:  %.2f GB" % (8.0 * nnz / 2**30))
    if assembly == "affine":
        print("  affine cache:   %.2f GB  (8 * (d+1) * nnz, d = %d)" % (8.0 * (d + 1) * nnz / 2**30, d))
    else:
        # measured: 1.65 s of FEniCS assembly at 160^2, essentially linear in the cell count
        print("  per-sample assembly: ~%.0f s/sample projected (1.0 s at 160^2, linear in cells)."
              % (1.0 * (free / 25281.0)))
        print("  -> prefer --assembly affine unless the cache does not fit")
    # measured on this machine: CHOLMOD factor ~2.25 kB/dof at 160^2 and 640^2, growing ~log(dofs)
    fac = 2.25e3 * free * (np.log(max(free, 3)) / np.log(408321.0))
    print("  a direct factor would be roughly %.1f GB (peak ~%.1f GB): %s"
          % (fac / 2**30, 1.9 * fac / 2**30,
             "fits in 32 GB" if 1.9 * fac < 32 * 2**30 else "DOES NOT fit in 32 GB, use --solver cg-amg"))


def converged(sampler):
    """None if the solver is exact, else (iterations, reason) with reason > 0 meaning success."""
    solver = sampler.solver
    if getattr(solver, "exact", True):
        return None
    ksp = getattr(solver, "_ksp", None)
    if ksp is None:
        return None
    return ksp.getIterationNumber(), ksp.getConvergedReason()


def draw(args, d):
    """Z as gt_meta defines it: one draw, shared with the graphing side."""
    return gt_meta.draw(args.nb_samples, d, args.seed)


def _measure(args, assembly, Z, exact, rtol, rebuild_every):
    """One row of the calibration: cg-amg at this tolerance and this reuse
    setting, against the exact values."""
    from SPDESampler.Solvers.linear_solvers import make_solver
    s = make_solver("cg-amg", rtol=rtol, atol=args.atol,
                    pc_rebuild_every=rebuild_every)
    sampler, _, _ = build(args, s, assembly)
    nit, bad, vals = [], 0, []
    t = time.time()
    for z in Z:
        vals.append(sampler.sample(z))
        c = converged(sampler)
        if c is not None:
            nit.append(c[0])
            if c[1] <= 0:
                bad += 1
    sec = (time.time() - t) / len(Z)
    del sampler
    dg = np.abs(np.array(vals) - exact)
    return {"sec": sec, "dmax": float(dg.max()),
            "rel": float((dg / np.abs(exact)).max()),
            "its": "%d-%d" % (min(nit), max(nit)) if nit else "-",
            "bad": bad}

# A banner for the calibration report alone. 
CAL_BANNER = """\
================================================================================
CALIBRATION OUTPUT -- there are NO reference values in this folder.

Produced by  compute_ground_truth.py --calibrate  (stage 1 of testing/launchGT.sh).
Its only purpose is to choose two settings for the real run: the solver tolerance
RTOL, and PC_REBUILD (how often the AMG hierarchy is rebuilt). Nothing downstream
reads this file -- it is for a human, and it is safe to delete at any time.

The actual ground truth lives in  GTresults_d<d>H<mesh>_<hash>/  -- never here.

  date      : %s
  identity  : %s
  mesh      : %d^2   <- NOT the mesh of the real run; a direct solve must fit here
  samples   : %d per row
================================================================================
"""

def calibrate(args):
    """Two questions that cannot be answered from theory, both measured on a mesh
    where an exact answer is still affordable:

      1. which rtol? The QoI error is bounded by ||w|| ||r||, not by rtol alone.
      2. can the AMG hierarchy be reused between samples, and how far?

    Both are answered against a direct solve on the SAME points, so the numbers
    in the 'max |dG|' column are the actual error of the reference, not a residual.

    TODO: Review this function. Append/correct if necessary, approve if ok.
    This has been fully vibe coded with Claude code and not reviewed yet.
    """
    n = args.calibrate_samples
    lines = []
    report = None
    if args.calibrate_out:
        ## Deliberately NOT a GTresults folder and deliberately not hash-named as
        ## one: these are advisory numbers, and the banner says so in case this
        ## folder is still lying around in six months. The hash is in the file
        ## name only so two coefficient sets cannot overwrite each other.
        idy = gt_meta.identity_from_args(args, args.mesh)
        os.makedirs(args.calibrate_out, exist_ok=True)
        report = os.path.join(args.calibrate_out,
                              "calibration_%s.txt" % gt_meta.gt_hash(idy))
        banner = CAL_BANNER % (time.strftime("%Y-%m-%d"),
                               gt_meta.canonical(idy), args.mesh, n)

    def emit(fmt="", *a):
        """Print AND rewrite the report after every line: this job can take hours
        and can hit its time limit, so the rows already measured must survive."""
        txt = (fmt % a) if a else fmt
        print(txt, flush=True)
        lines.append(txt)
        if report:
            with open(report, "w") as fh:
                fh.write(banner)
                fh.write("\n".join(lines) + "\n")

    emit("Calibration on a %d^2 mesh, %d samples per row" % (args.mesh, n))
    assembly = choose_assembly(args, 2 * args.nb_cosines)
    ref, disc, d = build(args, "cholmod", assembly)
    report_size(disc, d, assembly)
    Z = gt_meta.draw(n, d, args.seed)
    t = time.time()
    exact = np.array([ref.sample(z) for z in Z])
    direct = (time.time() - t) / n
    emit("  cholmod: %.2f s/solve, QoI mean %.12e, std %.3e\n"
         % (direct, exact.mean(), exact.std()))
    del ref

    emit("1) tolerance, rebuilding the preconditioner every sample")
    emit("   %-10s %-12s %-12s %-10s %-9s %s"
         % ("rtol", "max |dG|", "max |dG|/|G|", "s/solve", "vs direct", "CG its"))
    for rtol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-14):
        r = _measure(args, assembly, Z, exact, rtol, 1)
        emit("   %-10.0e %-12.3e %-12.3e %-10.3f %-9s %s%s"
             % (rtol, r["dmax"], r["rel"], r["sec"], "%.1fx" % (direct / r["sec"]),
                r["its"], "  DIVERGED x%d" % r["bad"] if r["bad"] else ""))
    emit("   -> pick the loosest rtol whose max |dG| sits at least two orders of magnitude")
    emit("      below the smallest quantity you report against this reference")
    emit("      (your tol_res, down to 5.7e-10). Put it in --rtol / RTOL.\n")

    emit("2) preconditioner reuse at rtol = %.0e (--pc_rebuild_every / PC_REBUILD)" % args.rtol)
    emit("   The AMG hierarchy is the expensive part and the sparsity pattern never")
    emit("   changes between samples, so reuse is the main speed lever. What it costs")
    emit("   is CG iterations -- and, if pushed too far, the tolerance itself.")
    emit("   %-14s %-12s %-12s %-10s %-9s %s"
         % ("rebuild every", "max |dG|", "max |dG|/|G|", "s/solve", "speedup", "CG its"))
    base = None
    for every in (1, 2, 5, 10, 25, 0):
        r = _measure(args, assembly, Z, exact, args.rtol, every)
        base = base or r["sec"]
        emit("   %-14s %-12.3e %-12.3e %-10.3f %-9s %s%s"
             % ("every sample" if every == 1 else ("once" if every == 0 else every),
                r["dmax"], r["rel"], r["sec"], "%.2fx" % (base / r["sec"]),
                r["its"], "  DIVERGED x%d" % r["bad"] if r["bad"] else ""))
    emit("   -> take the largest k whose max |dG| still matches row 1 at the same rtol")
    emit("      and that never DIVERGED. If no row is both faster and as accurate, leave")
    emit("      PC_REBUILD at 1: the reuse is not worth it at this mesh.")
    emit("")
    emit("These numbers are for %d^2. Costs scale, accuracy does not: the QoI error of"
         % args.mesh)
    emit("a given rtol is a property of the conditioning, which grows like h^-2, so a")
    emit("tolerance that is comfortable here is ~(mesh/%d)^2 times less so on the real" % args.mesh)
    emit("reference mesh. That is the argument for STRATEGY=rich.")

    if report:
        print("\nwrote %s" % report)


def yGT_file(args, out):
    return os.path.join(out, args.prefix + "yGT.npy")


def shard_file(args, out, shard=None):
    return os.path.join(out, args.prefix + "yGT.shard%dof%d.npy"
                        % (args.shard if shard is None else shard, args.nshards))


def drop_shards(out, start, end):
    n = 0
    for f in sorted(os.listdir(out)):
        if f.startswith(start) and f.endswith(end):
            os.remove(os.path.join(out, f))
            n += 1
    return n


def load_yGT(args, out, nb):
    """The authoritative values, as an array of length nb, NaN where missing.

    This file, and only this file, decides what is left to do. A shard must NOT
    look at the other shards' files when it works out its slice: on a cluster
    the tasks start minutes apart, so a later shard would see the earlier ones'
    progress, re-slice a smaller set, and leave gaps nobody computes.
    """
    vals = np.full(nb, np.nan)
    yf = yGT_file(args, out)
    if os.path.exists(yf):
        prev = np.load(yf)
        k = min(prev.size, nb)
        vals[:k] = prev[:k]
    return vals


def overlay_shards(args, out, vals, only=None):
    """Add the shard files' values to `vals`, in place; `only` restricts it to
    one file. They are scratch: each value sits next to the sample INDEX it
    belongs to, so deleting them costs re-solving and nothing else, and a round
    of 8 shards can pick up where a round of 16 stopped. Only finite values are
    overlaid -- a shard can add values, never remove them."""
    n = 0
    for f in sorted(os.listdir(out)):
        if not (f.startswith(args.prefix + "yGT.shard") and f.endswith(".npy")):
            continue
        if only is not None and os.path.join(out, f) != only:
            continue
        pack = np.load(os.path.join(out, f))
        if pack.ndim != 2:                 # not an (index, value) pack: ignore it
            continue
        i, v = pack[0].astype(int), pack[1]
        keep = (i < vals.size) & np.isfinite(v)
        vals[i[keep]] = v[keep]
        n += 1
    return n


def save_yGT(args, out, vals):
    """Write the reference, NEVER shrinking it.

    A merge against half the shards, or a run asking for fewer samples than are
    already there, used to overwrite this file with the leading complete block
    -- turning 2000 good values into a handful. The file is always nb_samples
    long now, with NaN where a sample is still missing, and the reader checks
    for those instead of inferring completeness from the length.
    """
    yf = yGT_file(args, out)
    if os.path.exists(yf):
        old = np.load(yf)
        if old.size > vals.size:
            merged = old.copy()
            new = np.isfinite(vals)
            merged[:vals.size][new] = vals[new]
            vals = merged
    np.save(yf, vals)
    return yf


def report_missing(args, vals, what):
    miss = np.flatnonzero(~np.isfinite(vals))
    if miss.size == 0:
        return
    print("%d of %d samples are still missing (%s): %s%s"
          % (miss.size, vals.size, what, miss[:20].tolist(),
             " ..." if miss.size > 20 else ""))
    print("Rerun to compute exactly these; loosen --rtol or use a direct "
          "--solver first if they are non-convergences.")
    if not args.allow_failures:
        sys.exit(1)


def merge(args, d, out):
    """Fold every shard file into the reference. Idempotent, and it cannot lose
    values: a missing shard leaves its samples NaN instead of truncating."""
    vals = load_yGT(args, out, args.nb_samples)
    n = overlay_shards(args, out, vals)
    yf = save_yGT(args, out, vals)
    print("folded in %d shard file(s)" % n)
    print("merged -> %s (%d of %d values)"
          % (yf, int(np.isfinite(vals).sum()), args.nb_samples))
    ## Every value the shard files held is now in yGT.npy, and that file is
    ## never shrunk, so they are redundant -- complete set or not. A shard that
    ## reruns derives its work from yGT.npy, not from its own file, so nothing
    ## is lost. Keeping them would also let a value computed at a looser --rtol
    ## be re-injected after a deliberate recompute with a direct solver, since
    ## the solver is (deliberately) not part of the identity.
    if n and not args.keep_shards:
        print("removed %d redundant shard file(s)"
              % drop_shards(out, args.prefix + "yGT.shard", ".npy"))
    report_missing(args, vals, "no shard file holds them")


def main(argv=None):
    args = apply_cfg(parse_args(argv), argv)
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")          # hypre and CHOLMOD are serial here; parallelise over samples

    if args.calibrate:
        return calibrate(args)

    d = 2 * args.nb_cosines
    ## The folder name carries a hash of the identity, so a set computed for a
    ## different coefficient (varsigma above all) cannot be resumed into or
    ## merged with this one. See utilities/gt_meta.py.
    identity = gt_meta.identity_from_args(args, args.mesh)
    out = args.out or gt_meta.folder(identity)
    os.makedirs(out, exist_ok=True)
    gt_meta.check(out, identity)               # exits if the folder is another problem
    print("identity: %s" % gt_meta.canonical(identity))
    print("folder:   %s" % out)

    def save_meta(n_done, its, failed):
        gt_meta.write(out, identity, {
            "solver": args.solver, "rtol": args.rtol, "atol": args.atol,
            "nb_samples_done": int(n_done),
            "max_iterations": (max(its) if its else None),
            "all_converged": not failed,
            "failed_samples": sorted(failed)})

    if args.merge:
        return merge(args, d, out)

    Z = draw(args, d)
    zf = os.path.join(out, args.prefix + "samples.npy")
    if os.path.exists(zf):
        ## Compare the OVERLAP, in both directions: the draw is prefix-stable,
        ## so growing a set must extend the stored points, never replace them.
        old = np.load(zf)
        k = min(old.shape[0], args.nb_samples)
        if old.shape[1] != Z.shape[1] or not np.allclose(old[:k], Z[:k]):
            sys.exit("%s holds different sample points: reuse its --seed or move it aside" % zf)
        if old.shape[0] > args.nb_samples:
            Z = old                       # never shrink the stored points either
    np.save(zf, Z)
    print("sample points: %s %s" % (zf, Z.shape))

    ## The work set is whatever <prefix>yGT.npy does not have yet: the tail,
    ## when nb_samples has grown, plus any NaN a failed solve left behind.
    ## Nothing depends on the shard files surviving -- see load_done.
    vals = load_yGT(args, out, args.nb_samples)
    todo_all = np.flatnonzero(~np.isfinite(vals))
    print("found %d of %d samples in %s"
          % (args.nb_samples - todo_all.size, args.nb_samples,
             os.path.basename(yGT_file(args, out))))

    if args.shard is not None:
        if args.nshards is None:
            sys.exit("--shard requires --nshards")
        ## The MISSING samples are what gets shared out, so a second round only
        ## splits the remainder. The slice is taken before any shard file is
        ## read, so every task of a round agrees on the same partition however
        ## far apart they start.
        yf = shard_file(args, out)
        mine = todo_all[args.shard::args.nshards]
        print("shard %d of %d: %d of the %d missing samples"
              % (args.shard, args.nshards, mine.size, todo_all.size))
        ## ...and only now my own file, to resume after a time limit
        overlay_shards(args, out, vals, only=yf)
        mine = np.array([k for k in mine if not np.isfinite(vals[k])], dtype=int)
        if mine.size:
            print("  %d of them are still to do" % mine.size)
    else:
        ## Unsharded: pick up a half-finished sharded round as well
        overlay_shards(args, out, vals)
        mine = np.flatnonzero(~np.isfinite(vals))
        yf = yGT_file(args, out)
    if mine.size == 0:
        print("nothing to do")
        return

    assembly = choose_assembly(args, d)
    sampler, disc, d = build(args, make(args), assembly)
    print("solver %s (rtol %.0e), %s assembly" % (args.solver, args.rtol, assembly))
    report_size(disc, d, assembly)

    failed = []

    def flush(n_attempted):
        seen = mine[:n_attempted]
        if args.shard is not None:
            ## (index, value) pairs, so this file stays readable whatever the
            ## work set looks like next time
            np.save(yf, np.vstack((seen, vals[seen])))
        else:
            save_yGT(args, out, vals)
        save_meta(int(np.isfinite(vals).sum()), its, failed)

    t0, its = time.time(), []
    for n, k in enumerate(mine, 1):
        try:
            vals[k] = sampler.sample(Z[k])
            c = converged(sampler)
            if c is not None:
                its.append(c[0])
                if c[1] <= 0:
                    raise RuntimeError("PETSc converged reason %d" % c[1])
        except Exception as exc:                            # noqa: BLE001
            vals[k] = np.nan
            failed.append(int(k))
            print("  !! sample %d failed: %s" % (k, exc), flush=True)
        if n % args.checkpoint == 0 or n == mine.size:
            flush(n)
            el = time.time() - t0
            print("  %d/%d  %.2f s/solve  eta %.1f min%s%s"
                  % (n, mine.size, el / n, (mine.size - n) * el / n / 60.0,
                     "  CG its %d-%d" % (min(its), max(its)) if its else "",
                     "  FAILED %d" % len(failed) if failed else ""), flush=True)
    flush(mine.size)
    print("wrote %s (%.1f min)" % (yf, (time.time() - t0) / 60.0))

    if args.shard is not None:
        print("run with --merge once every shard is done")
        report_missing(args, vals[mine], "this shard could not solve them")
        return
    report_missing(args, vals, "they did not converge")


if __name__ == "__main__":
    main()
