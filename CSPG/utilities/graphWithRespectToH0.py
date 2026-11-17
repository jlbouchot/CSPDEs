"""Compute time and error as the base mesh h0 = 1/mesh_x moves, at fixed L-J.

The sweep holds the TARGET mesh and the NUMBER OF LEVELS fixed, so every run
discretises on the same grids and recovers to the same tolerance; the only thing
that differs is the first-level sparsity the schedule asks for. The figures
therefore price the schedule, not the discretisation -- which is what makes this
the calibration counterpart to Exp1, where the grids and the schedule move
together and cannot be told apart.

Because the grids are shared, the invariants are checkable, and they are checked
rather than assumed: a row whose target mesh or level count drifted is not
comparable to the others and is refused here instead of being plotted.
"""

import os.path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import Check_ML
from utilities.utils import getComputeTimes, ground_truth_for, load_runs, savefig

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/11/17"
__created__ = "2026/08/14"


# ============================== settings ====================================
HERE = os.path.dirname(os.path.abspath(__file__))
TESTING = os.path.join(HERE, '..', 'testing')

# Meshes used for the trials. 
MESH_X = [5, 10, 20, 40, 80]
# Location of the results
RESULTS = os.path.join(TESTING, 'results', 'avgqoi_Exp4meshxvaries')

NB_TESTS = 2000

# Configuration for the ground truth computations.
GT_CFG = os.path.join(HERE, '..', 'data', 'Exp4_meshxvaries.ini')
GT_SEED = 20261106
GT_MESH = 'rich640-1280-2560'
# ============================================================================


# ------------------------------------------------------------------ load ---
if not os.path.isdir(RESULTS):
    raise SystemExit(
        "no Exp4 results in %s.\n"
        "Generate them first:\n"
        "    locally   cd testing && bash Exp4NoBatch.sh\n"
        "    on margaret  cd testing && bash launchExp4.sh\n"
        "The sweep needs no new ground truth: mesh_x is not part of the\n"
        "reference identity, so it scores against the one Exp1/Exp2/Exp3 use."
        % RESULTS)

runs = load_runs(RESULTS, dict((mx, 'mx_val_%d' % mx) for mx in MESH_X))

d = runs[MESH_X[0]].cspde_result[0].d

# -------------------------------------------------------- the invariants ---
shape = {}
for mx in MESH_X:
    run = runs[mx]
    stored_mx = run.spde_model[5][0]
    L = run.L
    nlev = len(run.cspde_result)
    J = L - nlev + 1
    if int(stored_mx) != int(mx):
        raise SystemExit(
            "mx_val_%d in %s was run with mesh_x = %s.\n"
            "The file name and the stored model disagree, so the sweep cannot "
            "be labelled." % (mx, RESULTS, stored_mx))
    shape[mx] = dict(mesh_x=int(stored_mx), J=J, L=L, nlev=nlev,
                     coarsest=int(stored_mx) * 2 ** J,
                     finest=int(stored_mx) * 2 ** L)

targets = set(s['finest'] for s in shape.values())
nlevels = set(s['nlev'] for s in shape.values())
if len(targets) != 1 or len(nlevels) != 1:
    for mx in MESH_X:
        s = shape[mx]
        print("  mesh_x=%-4d J=%d L=%d  %d levels  grids %d^2..%d^2"
              % (mx, s['J'], s['L'], s['nlev'], s['coarsest'], s['finest']))
    raise SystemExit(
        "the rows do not share a target mesh (%s) and a level count (%s).\n"
        "Exp4 prices the sparsity schedule by holding both fixed; without that "
        "the rows differ in tol_res and in their detail sparsities too, and the "
        "comparison means nothing. Fix the (mesh_x, l_start, nb_level) grid in "
        "testing/launchExp4.sh -- its own guards refuse this before submitting."
        % (sorted(targets), sorted(nlevels)))

TARGET = targets.pop()
NLEV = nlevels.pop()
print("\n  all %d rows: target %d^2, %d levels, grids %s"
      % (len(MESH_X), TARGET, NLEV,
         ", ".join("%d^2" % (shape[MESH_X[0]]['coarsest'] * 2 ** i)
                   for i in range(NLEV))))

# -------------------------------------------------------- ground truth -----
y_GT, Z, identity = ground_truth_for(runs, GT_CFG, GT_MESH, GT_SEED, d,
                                     NB_TESTS, base_dir=TESTING)

# ------------------------------------------------- errors and timings ------
linferror, rmserror, relinf = {}, {}, {}
t_total, t_pde, t_recovery = {}, {}, {}
cpu_total = {}
s_first, m_first = {}, {}

# Expected user/wall for a healthy run: numba at 4 threads measured 3.63x
# (see testing/launchExp3.sh). A row that ran while the node was busy shows a
# MUCH lower ratio
CPU_PER_WALL_EXPECTED = 3.63
CPU_PER_WALL_FLOOR = 2.5

for mx in MESH_X:
    run = runs[mx]
    err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
    linferror[mx] = np.abs(err).max()
    rmserror[mx] = np.linalg.norm(err) / np.sqrt(err.size)
    relinf[mx] = linferror[mx] / np.abs(y_GT).max()

    s_first[mx] = run.cspde_result[0].s
    m_first[mx] = run.cspde_result[0].m

    pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
    t_total[mx] = pde[0] + recovery[0] + mtx[0] + js[0]
    t_pde[mx] = pde[0]
    t_recovery[mx] = recovery[0]
    cpu_total[mx] = (pde[1] + recovery[1] + mtx[1] + js[1]
                     + pde[2] + recovery[2] + mtx[2] + js[2])

print("\n  mesh_x  h0      J  L   s_J    m_J       linf        rel linf    rms         "
      "wall (h)   CPU (h)   cpu/wall  PDE (h)   recovery (h)  rec share")
for mx in MESH_X:
    s = shape[mx]
    share = t_recovery[mx] / t_total[mx] if t_total[mx] > 0 else float('nan')
    ratio = cpu_total[mx] / t_total[mx] if t_total[mx] > 0 else float('nan')
    print("  %-7d 1/%-5d %-2d %-3d %-6g %-9g %-11.4g %-11.4g %-11.4g %-10.3f %-9.3f %-9.2f %-9.3f %-13.3f %.0f%%"
          % (mx, mx, s['J'], s['L'], s_first[mx], m_first[mx],
             linferror[mx], relinf[mx], rmserror[mx],
             t_total[mx] / 3600, cpu_total[mx] / 3600, ratio,
             t_pde[mx] / 3600, t_recovery[mx] / 3600, 100 * share))
print()

## ------------------------------------------------------ contention check ---
busy = [mx for mx in MESH_X
        if t_total[mx] > 0 and cpu_total[mx] / t_total[mx] < CPU_PER_WALL_FLOOR]
if busy:
    print("  WARNING: cpu/wall below %.1f (expected ~%.2f at numba=4) for mesh_x in %s."
          % (CPU_PER_WALL_FLOOR, CPU_PER_WALL_EXPECTED, busy))
    print("           Those rows spent wall clock WAITING, not computing: they "
          "shared a node.\n           Compare the CPU (h) column across rows "
          "instead of wall (h), or re-run\n           those rows on a quiet node "
          "before quoting any wall-clock ratio.")
    print()

# ------------------------------------------------------------- null test ---
by_plan = {}
for mx in MESH_X:
    by_plan.setdefault(s_first[mx], []).append(mx)
for s_val, group in sorted(by_plan.items()):
    if len(group) < 2:
        continue
    errs = [linferror[mx] for mx in group]
    # CPU-seconds, not wall clock: two identical rows scheduled against
    # different node loads differ in wall clock for reasons that have nothing
    # to do with the method, which would make a sound null test look broken.
    cpus = [cpu_total[mx] for mx in group]
    spread_e = (max(errs) - min(errs)) / max(min(errs), 1e-300)
    spread_t = (max(cpus) - min(cpus)) / max(min(cpus), 1e-300)
    verdict = "OK" if spread_e < 1e-9 else ("close" if spread_e < 0.05 else "!! DIFFER !!")
    print("  null test, s_J = %g, mesh_x in %s:  error spread %.2e, "
          "CPU-time spread %.1f%%  -> %s"
          % (s_val, group, spread_e, 100 * spread_t, verdict))
    if verdict == "!! DIFFER !!":
        print("      Identical plans on identical grids gave different errors. "
              "Something outside\n      CSPDE_v2.level_plan depends on the level "
              "index; do not read the trend below.")
print()

# ----------------------------------------------------------- figures -------
h0 = [1.0 / mx for mx in MESH_X]
colours = plt.cm.rainbow(np.linspace(0, 1, 4))
markers = ["o", "v", "^", "<"]

## 1. the headline: cost and error against h0, with the grids held fixed.
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'$h_0 = 1/\mathrm{mesh}_x$   '
               '(target $%d^2$, %d levels, identical grids)' % (TARGET, NLEV))
ax1.set_ylabel('Wall clock time (h)')
for i, (series, label) in enumerate(((t_total, 'Total time'),
                                     (t_pde, 'PDE solve time'),
                                     (t_recovery, 'Sparse recovery time'))):
    ax1.plot(h0, [series[mx] / 3600 for mx in MESH_X],
             color=colours[i], marker=markers[i], label=label)
ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xticks(h0)
ax1.set_xticklabels(['1/%d' % mx for mx in MESH_X])

ax2 = ax1.twinx()
ax2.set_ylabel(r'$\ell_\infty$ error')
ax2.set_yscale('log')
ax2.plot(h0, [linferror[mx] for mx in MESH_X],
         color=colours[3], marker=markers[3], label='Final approximation error')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)
fig.tight_layout()
savefig("timesAndErrorWRTh0")

## 2. the same against what actually varies: the first-level sparsity. The
## x-axis here is the schedule's own ask, so a flat error curve is a direct
## statement that the const_sj floor over-provisions.
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'first-level sparsity $s_J$ asked for by the schedule')
ax1.set_ylabel('Wall clock time (h)')
order = sorted(MESH_X, key=lambda mx: s_first[mx])
ax1.plot([s_first[mx] for mx in order], [t_total[mx] / 3600 for mx in order],
         color=colours[0], marker='o', label='Total time')
ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax2 = ax1.twinx()
ax2.set_ylabel(r'$\ell_\infty$ error')
ax2.set_yscale('log')
ax2.plot([s_first[mx] for mx in order], [linferror[mx] for mx in order],
         color=colours[3], marker='<', label='Final approximation error')
for mx in order:
    ax2.annotate(r'$h_0=1/%d$' % mx, (s_first[mx], linferror[mx]),
                 textcoords='offset points', xytext=(4, 4), fontsize=8)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)
fig.tight_layout()
savefig("errorAndTimeWRTsJ")

## 3. where the time goes. The recovery grows with s_J while the PDE work is
## nearly fixed (the grids never change), so this is the shape of the price.
plt.figure()
idx = np.arange(len(MESH_X))
pde_share = np.array([t_pde[mx] / t_total[mx] for mx in MESH_X])
rec_share = np.array([t_recovery[mx] / t_total[mx] for mx in MESH_X])
plt.bar(idx, pde_share, color=colours[1], label='PDE solves')
plt.bar(idx, rec_share, bottom=pde_share, color=colours[2], label='Sparse recovery')
plt.bar(idx, 1.0 - pde_share - rec_share, bottom=pde_share + rec_share,
        color=colours[0], label='Matrix + index set')
plt.xticks(idx, [r'$1/%d$' % mx + '\n$s_J=%g$' % s_first[mx] for mx in MESH_X])
plt.xlabel(r'$h_0$')
plt.ylabel('share of the total wall clock')
plt.ylim(0, 1)
plt.legend()
plt.tight_layout()
savefig("timeSharesWRTh0")

## 4. work-precision. Every point buys its accuracy on the SAME grids, so a
## horizontal spread here is pure schedule cost.
##
## The x-axis is CPU-seconds, not wall clock, ON PURPOSE. Rows submitted at
## different times meet different node loads, and wall clock absorbs that
## waiting; CPU time does not. Since the horizontal gap IS the result here,
## it must not be a scheduling artefact. Swap to t_total if you ever need the
## wall-clock version for a cost-to-the-user argument, but check the cpu/wall
## column first.
plt.figure()
handles = []
for i, mx in enumerate(MESH_X):
    handles.append(plt.scatter(cpu_total[mx], linferror[mx],
                               color=plt.cm.rainbow(i / max(len(MESH_X) - 1, 1)),
                               marker=markers[i % len(markers)]))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('CPU time, user + system (s)')
plt.ylabel(r'$\ell_\infty$ error')
plt.legend(handles=handles,
           labels=[r'$h_0 = 1/%d$  ($s_J = %g$)' % (mx, s_first[mx]) for mx in MESH_X])
plt.tight_layout()
savefig("errorVsWorkWRTh0")

plt.show()
