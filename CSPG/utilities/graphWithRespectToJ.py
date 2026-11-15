"""Exp3: compute time and error as the coarsest level J moves, with L fixed.

Reads one experiment directory of ``J_val_<J>`` shelves. The target level L is
the same for every point, so all of them share ONE reference and differ only in
how much of the work is pushed down onto coarse meshes. J = L is the
single-level approach, which is why that point doubles as the ML-vs-SL
comparison.

The QoI the runs actually computed is read back from the shelves and checked
against the reference before anything is plotted: scoring against the wrong
functional yields a smooth, plausible curve rather than an error.
"""

import os.path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import Check_ML                 # noqa: F401 -- defines the namedtuples the shelves unpickle into
from utilities.utils import getComputeTimes, ground_truth_for, load_runs

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/11/14"


# ============================== settings ====================================
HERE = os.path.dirname(os.path.abspath(__file__))
TESTING = os.path.join(HERE, '..', 'testing')

## The sweep. L is FIXED (testing/launchExp3.sh: L_VAL), so J = L_VAL is the
## single-level reference and small J means more detail levels. Below J = 3 the
## first-level sparsity stops being flat and the runs become infeasible, which
## is why the grid stops there.
JS = [2, 3, 4, 5, 6]
L_VAL = 6
## See this folder's README.QOI-MISHAP.md: these runs computed 'avg'.
RESULTS = os.path.join(TESTING, 'results', 'avgqoi_Exp3Jvaries_preQoIfix')

NB_TESTS = 2000

## Configuration for the ground truth.
GT_CFG = os.path.join(HERE, '..', 'data', 'Exp3_Jvaries.ini')
GT_SEED = 20261106
# Mesh type used for computing the ground truth.
GT_MESH = 'rich640-1280-2560'
# ============================================================================


# ------------------------------------------------------------------ load ---
runs = load_runs(RESULTS, dict((J, 'J_val_%d' % J) for J in JS))

d = runs[JS[0]].cspde_result[0].d
base_mesh = runs[JS[0]].spde_model[5][0] \
    if isinstance(runs[JS[0]].spde_model, (list, tuple)) else None

# -------------------------------------------------------- ground truth -----
y_GT, Z, identity = ground_truth_for(runs, GT_CFG, GT_MESH, GT_SEED, d,
                                     NB_TESTS, base_dir=TESTING)

# ------------------------------------------------- errors and timings ------
linferror, rmserror, relinf = {}, {}, {}
t_total, t_pde, t_recovery = {}, {}, {}
t_user, t_sys = {}, {}

for J in JS:
    run = runs[J]
    err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
    linferror[J] = np.abs(err).max()
    rmserror[J] = np.linalg.norm(err) / np.sqrt(err.size)
    relinf[J] = linferror[J] / np.abs(y_GT).max()

    pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
    ## Each is [wall clock, user, system]; seconds throughout.
    t_total[J] = pde[0] + recovery[0] + mtx[0] + js[0]
    t_pde[J] = pde[0]
    t_recovery[J] = recovery[0]
    t_user[J] = pde[1] + recovery[1] + mtx[1] + js[1]
    t_sys[J] = pde[2] + recovery[2] + mtx[2] + js[2]

best = min(JS, key=lambda J: t_total[J])
print("\n  J   levels  coarsest  linf        rel linf    rms         "
      "total (h)  PDE (h)   recovery (h)  user (h)  sys (h)")
for J in JS:
    coarsest = '%d^2' % (base_mesh * 2 ** J) if base_mesh else '?'
    print("  %-3d %-7d %-9s %-11.4g %-11.4g %-11.4g %-10.3f %-9.3f %-13.3f %-9.3f %.3f%s"
          % (J, L_VAL - J + 1, coarsest, linferror[J], relinf[J], rmserror[J],
             t_total[J] / 3600, t_pde[J] / 3600, t_recovery[J] / 3600,
             t_user[J] / 3600, t_sys[J] / 3600,
             "   <- cheapest" if J == best else ""))
print("\n  single level (J = L = %d): %.3f h;  cheapest (J = %d): %.3f h"
      % (L_VAL, t_total[L_VAL] / 3600, best, t_total[best] / 3600), end="")
if t_total[best] > 0:
    print("  -> %.1fx" % (t_total[L_VAL] / t_total[best]))
else:
    print()
print()

# ----------------------------------------------------------- figures -------
colours = plt.cm.rainbow(np.linspace(0, 1, 4))
markers = ["o", "v", "^", "<"]

fig, ax1 = plt.subplots()
ax1.set_xlabel('J  (L = %d; J = %d is single level)' % (L_VAL, L_VAL))
ax1.set_ylabel('Wall clock time (h)')
for i, (series, label) in enumerate(((t_total, 'Total time'),
                                     (t_pde, 'PDE solve time'),
                                     (t_recovery, 'Sparse recovery time'))):
    ax1.plot(JS, [series[J] / 3600 for J in JS],
             color=colours[i], marker=markers[i], label=label)
ax1.set_xticks(JS)

ax2 = ax1.twinx()
ax2.set_ylabel(r'$\ell_\infty$ error')
ax2.set_yscale('log')
ax2.plot(JS, [linferror[J] for J in JS],
         color=colours[3], marker=markers[3], label='Final approximation error')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)
fig.tight_layout()
plt.savefig("timesAndErrorWRTjs.eps")
plt.savefig("timesAndErrorWRTjs.jpg")

if base_mesh:
    coarse = [base_mesh * 2 ** J for J in JS]
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Coarsest grid of the hierarchy (points per side; '
                   r'$h_J = 1/\mathrm{width}$, target $%d^2$)' % (base_mesh * 2 ** L_VAL))
    ax1.set_ylabel('Wall clock time (h)')
    for i, (series, label) in enumerate(((t_total, 'Total time'),
                                         (t_pde, 'PDE solve time'),
                                         (t_recovery, 'Sparse recovery time'))):
        ax1.plot(coarse, [series[J] / 3600 for J in JS],
                 color=colours[i], marker=markers[i], label=label)
    ax1.set_xscale('log', base=2)
    ax1.set_xticks(coarse)
    ax1.set_xticklabels([str(c) for c in coarse])

    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\ell_\infty$ error')
    ax2.set_yscale('log')
    ax2.plot(coarse, [linferror[J] for J in JS],
             color=colours[3], marker=markers[3], label='Final approximation error')
    for J in JS:
        ax2.annotate('J=%d (%d lvl)' % (J, L_VAL - J + 1),
                     (base_mesh * 2 ** J, linferror[J]),
                     textcoords='offset points', xytext=(4, 4), fontsize=8)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2)
    fig.tight_layout()
    plt.savefig("timesAndErrorWRTcoarsestMesh.eps")
    plt.savefig("timesAndErrorWRTcoarsestMesh.jpg")

plt.figure()
handles = []
for i, J in enumerate(JS):
    handles.append(plt.scatter(t_total[J], linferror[J],
                               color=plt.cm.rainbow(i / max(len(JS) - 1, 1)),
                               marker=markers[i % len(markers)]))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wall clock computing time (s)')
plt.ylabel(r'$\ell_\infty$ error')
plt.legend(handles=handles,
           labels=["J = %d%s" % (J, " (single level)" if J == L_VAL else "")
                   for J in JS])
plt.tight_layout()
plt.savefig("errorVsWorkWRTjs.eps")
plt.savefig("errorVsWorkWRTjs.jpg")

plt.show()
