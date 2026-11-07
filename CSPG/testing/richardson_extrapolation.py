"""Build a reference QoI beyond the mesh that fits in memory, by extrapolating in h.

The QoI of this problem converges at order 2 in the mesh width: the measured detail magnitude
RMS(G(2n) - G(n)) falls by a factor 4.00 per refinement. So with

    G(n) = G* + C n^-q + O(n^-(q+1)),      q ~ 2,

two meshes already give G* to one extra order, and three let the order be measured rather than
assumed:

    q_hat  = log2( (G(n) - G(2n)) / (G(2n) - G(4n)) )
    G*     = G(4n) + ( G(4n) - G(2n) ) / (2^q_hat - 1).

Why this is worth doing: a direct solve on a 3000^2 mesh needs a ~25 GB Cholesky factor and still
leaves an O(n^-2) error. Extrapolating from 640/1280/2560 peaks at a 12 GB affine cache with
cg-amg, and leaves an O(n^-3) error -- the accuracy a single P1 solve would only reach at roughly
n = 2560^1.5, i.e. about 130000^2. The script reports the measured correction and a
self-consistency estimate instead of asking anyone to believe that.

It writes the same two files as compute_ground_truth.py, in a folder named by the same
utilities/gt_meta.py identity hash (with "rich640-1280-2560" in place of the single mesh width),
so the two kinds of reference cannot be mixed up:

    <out>/<prefix>samples.npy, <out>/<prefix>yGT.npy

and additionally keeps the per-mesh values in <out>/<prefix>levels.npz so the extrapolation can be
redone, or its order inspected, without re-solving.

Usage
-----
    python richardson_extrapolation.py --meshes 640 1280 2560 --nb_samples 1000 --solver cg-amg
    python richardson_extrapolation.py --meshes 640 1280 2560 --nb_samples 1000 --shard 3 --nshards 8
    python richardson_extrapolation.py --meshes 640 1280 2560 --nb_samples 1000 --merge
"""

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019 - 2026, INRIA, LMU Munisch, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/11/03"
__lastmodified__ = "2026/11/07"

import argparse
import copy
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np

from compute_ground_truth import (add_common_options, apply_cfg, build, choose_assembly,
                                  converged, draw, drop_shards, make, report_missing,
                                  report_size, save_yGT, yGT_file)
from utilities import gt_meta


def parse_args_parser():
    """The shared options come from compute_ground_truth.add_common_options, so
    the two scripts cannot drift apart; only the two options that are genuinely
    about extrapolation are declared here."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_options(p)
    p.add_argument("--meshes", type=int, nargs="+", default=[640, 1280, 2560],
                   help="two or more meshes, coarsest first, each twice the previous "
                        "(two give order-2 Richardson, three or more let the order be measured)")
    p.add_argument("--order", type=float, default=None,
                   help="force the convergence order instead of measuring it (2.0 is the theory)")
    return p


def parse_args(argv=None):
    return parse_args_parser().parse_args(argv)


def levels_file(args, out):
    return os.path.join(out, args.prefix + "levels.npz")


def levels_shard_file(args, out):
    return os.path.join(out, args.prefix + "levels.shard%dof%d.npz"
                        % (args.shard, args.nshards))


def load_levels(args, out, meshes, nb):
    """The authoritative per-mesh values, (nb_meshes, nb), NaN where missing.

    As in compute_ground_truth.load_yGT, a shard works out its slice from THIS
    file alone and never from the other shards' files: tasks start minutes
    apart, so a later shard that saw the earlier ones' progress would re-slice
    a smaller set and leave gaps nobody computes.
    """
    G = np.full((len(meshes), nb), np.nan)
    f = levels_file(args, out)
    if os.path.exists(f):
        z = np.load(f)
        if list(np.asarray(z["meshes"]).tolist()) != list(meshes):
            sys.exit("%s was computed on meshes %s, not %s"
                     % (f, np.asarray(z["meshes"]).tolist(), meshes))
        prev = z["G"]
        k = min(prev.shape[1], nb)
        G[:, :k] = prev[:, :k]
    return G


def overlay_level_shards(args, out, G, only=None):
    """Add the shard files' columns to G, in place; `only` restricts it to one
    file. Each column carries its sample index, so the files are scratch:
    deleting them costs re-solving and nothing else."""
    n = 0
    for f in sorted(os.listdir(out)):
        if not (f.startswith(args.prefix + "levels.shard") and f.endswith(".npz")):
            continue
        path = os.path.join(out, f)
        if only is not None and path != only:
            continue
        z = np.load(path)
        idx, g = z["idx"].astype(int), z["G"]
        if g.shape[0] != G.shape[0] or g.shape[1] != idx.size:
            continue
        keep = idx < G.shape[1]
        tgt, sub = idx[keep], g[:, keep]
        G[:, tgt] = np.where(np.isfinite(sub), sub, G[:, tgt])
        n += 1
    return n


def save_levels(args, out, G, meshes):
    """Write the per-mesh values, NEVER shrinking: a run for fewer samples must
    not drop columns a bigger run computed."""
    f = levels_file(args, out)
    if os.path.exists(f):
        old = np.load(f)["G"]
        if old.shape[1] > G.shape[1]:
            wide = old.copy()
            fin = np.isfinite(G)
            wide[:, :G.shape[1]] = np.where(fin, G, wide[:, :G.shape[1]])
            G = wide
    np.savez(f, G=G, meshes=np.array(meshes), idx=np.arange(G.shape[1]))
    return G


def extrapolate(G, order=None):
    """G has shape (nb_meshes, nb_samples), coarsest first. Returns (G_star, q_hat, diagnostics)."""
    G = np.asarray(G, dtype=float)
    ## The order and the error bar are measured on the COMPLETE samples only;
    ## G_star below keeps the full width, so an incomplete sample stays NaN.
    Gc = G[:, np.isfinite(G).all(axis=0)]
    d_fine = G[-1] - G[-2]
    if G.shape[0] >= 3 and order is None:
        d_coarse = Gc[-2] - Gc[-3]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = d_coarse / (Gc[-1] - Gc[-2])
        ok = np.isfinite(ratio) & (ratio > 0)
        # one order for the whole set: per-sample ratios are noisy where the detail nearly vanishes
        q = float(np.log2(np.median(ratio[ok]))) if ok.any() else 2.0
    else:
        q = 2.0 if order is None else float(order)
    G_star = G[-1] + d_fine / (2.0 ** q - 1.0)
    diag = {"order": q,
            "correction_rms": float(np.sqrt(np.mean(((Gc[-1] - Gc[-2]) / (2.0 ** q - 1.0)) ** 2))),
            "finest_rms": float(np.sqrt(np.mean(Gc[-1] ** 2)))}
    if G.shape[0] >= 3:
        # redo the extrapolation one level coarser: the gap between the two is an error estimate
        G_star = Gc[-1] + (Gc[-1] - Gc[-2]) / (2.0 ** q - 1.0)
        G_star_coarse = Gc[-2] + (Gc[-2] - Gc[-3]) / (2.0 ** q - 1.0)
        diag["selfcheck_rms"] = float(np.sqrt(np.mean((G_star - G_star_coarse) ** 2)))
        diag["selfcheck_max"] = float(np.max(np.abs(G_star - G_star_coarse)))
    ## full width again, so the caller's indices still line up
    G_star = G[-1] + d_fine / (2.0 ** q - 1.0)
    return G_star, q, diag


def report(diag, meshes):
    print("\nExtrapolation from meshes %s" % (meshes,))
    print("  measured order q           : %.4f   (theory: 2)" % diag["order"])
    print("  RMS of the finest QoI      : %.6e" % diag["finest_rms"])
    print("  RMS of the correction      : %.6e  (%.3g of the QoI)"
          % (diag["correction_rms"], diag["correction_rms"] / diag["finest_rms"]))
    if "selfcheck_rms" in diag:
        print("  RMS of the same extrapolation one level coarser differs by %.3e (max %.3e)"
              % (diag["selfcheck_rms"], diag["selfcheck_max"]))
        print("  -> take that as the reference's own error bar; it should sit well below any")
        print("     tolerance you report against this reference.")
    else:
        print("  pass three meshes to get an error bar on the extrapolation itself.")


def main(argv=None):
    args = apply_cfg(parse_args(argv), argv)
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")

    meshes = sorted(args.meshes)
    if len(meshes) < 2:
        sys.exit("--meshes needs at least two values")
    for a, b in zip(meshes, meshes[1:]):
        if b != 2 * a:
            print("WARNING: %d -> %d is not a factor-2 refinement; the extrapolation assumes it is"
                  % (a, b))

    d = 2 * args.nb_cosines
    ## Same identity, and therefore the same hash, as compute_ground_truth.py:
    ## an extrapolated reference is interchangeable with a single-mesh one only
    ## if the coefficient, QoI and seed agree. The "mesh" field spells out the
    ## ladder instead of one width (and the forced order, when there is one,
    ## since it changes the values) -- that keeps IDENTITY itself unchanged.
    tag = "rich" + "-".join(str(m) for m in meshes)
    if args.order is not None:
        tag += "q%g" % args.order
    identity = gt_meta.identity_from_args(args, tag)
    out = args.out or gt_meta.folder(identity)
    os.makedirs(out, exist_ok=True)
    gt_meta.check(out, identity)           # exits if the folder holds another problem
    print("identity: %s" % gt_meta.canonical(identity))
    print("folder:   %s" % out)
    lvl_file = os.path.join(out, args.prefix + "levels.npz")

    def save_meta(G, failed, diag=None):
        q = dict(solver=args.solver, rtol=args.rtol, atol=args.atol,
                 meshes=list(meshes),
                 nb_samples_done=int(np.isfinite(G).all(axis=0).sum()),
                 all_converged=not failed, failed_samples=sorted(set(failed)))
        if diag:
            q.update({"measured_order": diag["order"],
                      "correction_rms": diag["correction_rms"],
                      "selfcheck_rms": diag.get("selfcheck_rms")})
        gt_meta.write(out, identity, q)

    Z = draw(args, d)
    np.save(os.path.join(out, args.prefix + "samples.npy"), Z)

    def finish(G, failed, where):
        """Extrapolate what is complete, write both files, say what is missing."""
        G_star, q, diag = extrapolate(G, args.order)
        report(diag, meshes)
        save_meta(G, failed, diag)
        yf = save_yGT(args, out, G_star)
        print("\nwrote %s (%d of %d values)"
              % (yf, int(np.isfinite(G_star).sum()), args.nb_samples))
        report_missing(args, G_star, where)

    if args.merge:
        G = load_levels(args, out, meshes, args.nb_samples)
        n = overlay_level_shards(args, out, G)
        G = save_levels(args, out, G, meshes)
        print("folded in %d shard file(s)" % n)
        bad = np.flatnonzero(~np.isfinite(G).all(axis=0))
        ## redundant now that GT_levels.npz holds them; see merge() over there
        if n and not args.keep_shards:
            print("removed %d redundant shard file(s)"
                  % drop_shards(out, args.prefix + "levels.shard", ".npz"))
        finish(G[:, :args.nb_samples], [(int(i), 0) for i in bad],
               "no shard file holds them")
        return

    ## The work set is whatever GT_levels.npz is missing: the tail, when
    ## nb_samples has grown, plus any NaN a failed solve left behind. Nothing
    ## depends on the shard files surviving.
    G = load_levels(args, out, meshes, args.nb_samples)
    todo_all = np.flatnonzero(~np.isfinite(G).all(axis=0))
    print("found %d of %d samples in %s"
          % (args.nb_samples - todo_all.size, args.nb_samples,
             os.path.basename(levels_file(args, out))))

    if args.shard is not None:
        if args.nshards is None:
            sys.exit("--shard requires --nshards")
        ## The slice is taken before any shard file is read, so every task of a
        ## round agrees on the same partition however far apart they start.
        lvl_out = levels_shard_file(args, out)
        mine = todo_all[args.shard::args.nshards]
        print("shard %d of %d: %d of the %d missing samples"
              % (args.shard, args.nshards, mine.size, todo_all.size))
        overlay_level_shards(args, out, G, only=lvl_out)   # resume after a time limit
    else:
        lvl_out = levels_file(args, out)
        overlay_level_shards(args, out, G)                 # pick up a half-done round
        mine = np.flatnonzero(~np.isfinite(G).all(axis=0))
    if mine.size == 0:
        print("nothing to do")
        return

    def checkpoint():
        if args.shard is not None:
            np.savez(lvl_out, G=G[:, mine], meshes=np.array(meshes), idx=mine)
        else:
            save_levels(args, out, G, meshes)

    failed = []
    # one mesh at a time: the coarse levels are cheap and the fine level owns the memory,
    # so nothing is held open while the next discretisation is built
    for mi, mesh in enumerate(meshes):
        todo = np.array([k for k in mine if not np.isfinite(G[mi, k])], dtype=int)
        if todo.size == 0:
            print("mesh %d: already complete" % mesh)
            continue
        sub = copy.copy(args)
        sub.mesh = mesh
        assembly = choose_assembly(sub, d)
        sampler, disc, _ = build(sub, make(sub), assembly)
        print("\nmesh %d^2, %s assembly, solver %s: %d samples"
              % (mesh, assembly, args.solver, todo.size))
        report_size(disc, d, assembly)
        t0, its = time.time(), []
        for c, k in enumerate(todo):
            ## A sample that fails stays NaN and the run goes on: one bad solve
            ## must not throw away the hours already spent on the others, and
            ## NaN is also what the work set above looks for, so a later run
            ## retries exactly these.
            try:
                G[mi, k] = sampler.sample(Z[k])
                cv = converged(sampler)
                if cv is not None:
                    its.append(cv[0])
                    if cv[1] <= 0:
                        raise RuntimeError("PETSc converged reason %d" % cv[1])
            except Exception as exc:
                G[mi, k] = np.nan
                failed.append((int(k), mesh))
                print("  !! sample %d on mesh %d failed: %s" % (k, mesh, exc), flush=True)
            if (c + 1) % args.checkpoint == 0 or c + 1 == todo.size:
                checkpoint()
                el = time.time() - t0
                print("  %d/%d  %.2f s/solve  eta %.1f min%s"
                      % (c + 1, todo.size, el / (c + 1),
                         (todo.size - c - 1) * el / (c + 1) / 60.0,
                         "  CG its %d-%d" % (min(its), max(its)) if its else ""), flush=True)
        del sampler, disc

    checkpoint()
    if args.shard is not None:
        save_meta(G, failed)
        print("\nwrote %s; run with --merge once every shard is done" % lvl_out)
        report_missing(args, np.where(np.isfinite(G[:, mine]).all(axis=0), 0.0, np.nan),
                       "this shard could not solve them")
        return
    finish(G, failed, "they did not converge")


if __name__ == "__main__":
    main()
