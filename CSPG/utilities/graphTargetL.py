"""Exp1: approximation error and compute time against the target level L.

Reads one experiment directory of ``L_val_<L>`` shelves and plots them against
the ground truth that belongs to them.

A reference is valid for exactly one coefficient set, QoI, mesh and seed (see
``utilities/gt_meta.py``), and getting the QoI wrong is the dangerous case: it
produces a smooth, plausible-looking curve made of nothing but the gap between
two different functionals.  So the QoI is read back from the stored runs and
checked against the reference before anything is plotted.
"""

import os.path
import shelve
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import Check_ML
from utilities.utils import getComputeTimes, ground_truth_for, load_run, load_runs, savefig

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/11/16"


# ============================== settings ====================================
HERE = os.path.dirname(os.path.abspath(__file__))
TESTING = os.path.join(HERE, '..', 'testing')

# The runs to plot. Each is <RESULTS>/L_val_<L>; 
LMAX = [2, 3, 4, 5, 6, 7]
# Path to results
RESULTS = os.path.join(TESTING, 'results', 'avgqoi_h0_nbDisc_Fixed_preQoIfix')

# OPTIONAL single-level overlay, for the ML-vs-SL comparison.
#
# Simply adapt SL_RESULTS to point to a directory of shelves produced by the single-level sweeps.
# Make sure the qoi, coefficients, ... and all parameters match the experiments. 
# SL_RESULTS = None
SL_RESULTS = os.path.join(TESTING, 'results', 'avgqoi_SingleLevel')

NB_TESTS = 2000

# Configurations for the ground truth. 
GT_CFG = os.path.join(HERE, '..', 'data', 'Exp1_hfvaries.ini')
GT_SEED = 20261106
# mesh used for computing the reference.
GT_MESH = 'rich640-1280-2560'
# ============================================================================


# ------------------------------------------------------------------ load ---
results = load_runs(RESULTS, dict((L, 'L_val_%d' % L) for L in LMAX))

finest = results[LMAX[-1]]
d = finest.cspde_result[0].d
base_mesh = finest.spde_model[5][0] if isinstance(finest.spde_model, (list, tuple)) else None

# ------------------------------------------- optional single-level runs ----
sl_results = {}
if SL_RESULTS:
    if not os.path.isdir(SL_RESULTS):
        print("note: SL_RESULTS is set but %s does not exist; plotting the "
              "multilevel runs alone.\n      Generate them with: cd testing && "
              "bash launchSL.sh" % SL_RESULTS)
    else:
        for L in LMAX:
            if os.path.exists(os.path.join(SL_RESULTS, 'SL_val_%d.db' % L)) or \
               os.path.exists(os.path.join(SL_RESULTS, 'SL_val_%d.dat' % L)) or \
               os.path.exists(os.path.join(SL_RESULTS, 'SL_val_%d' % L)):
                sl_results[L] = load_run(SL_RESULTS, 'SL_val_%d' % L)
        if not sl_results:
            print("note: no SL_val_* shelves in %s; plotting the multilevel "
                  "runs alone." % SL_RESULTS)
        else:
            print("single level: L = %s"
                  % ", ".join(str(L) for L in sorted(sl_results)))

SL_LMAX = sorted(sl_results)

# -------------------------------------------------------- ground truth -----
both = dict(('ML L=%d' % L, results[L]) for L in LMAX)
both.update(('SL L=%d' % L, sl_results[L]) for L in SL_LMAX)
y_GT, Z, identity = ground_truth_for(both, GT_CFG, GT_MESH, GT_SEED, d,
                                     NB_TESTS, base_dir=TESTING)

# ------------------------------------------------- errors and timings ------
linferror, rmserror, relinf = {}, {}, {}
t_total, t_pde, t_recovery = {}, {}, {}
t_user, t_sys = {}, {}

for L in LMAX:
    run = results[L]
    err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
    linferror[L] = np.abs(err).max()
    rmserror[L] = np.linalg.norm(err) / np.sqrt(err.size)   # per sample, so it
    relinf[L] = linferror[L] / np.abs(y_GT).max()           # compares across NB_TESTS

    pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
    ## [0] wall clock, [1] user, [2] system. Seconds throughout; the plots
    ## convert where they need hours.
    t_total[L] = pde[0] + recovery[0] + mtx[0] + js[0]
    t_pde[L] = pde[0]
    t_recovery[L] = recovery[0]
    t_user[L] = pde[1] + recovery[1] + mtx[1] + js[1]
    t_sys[L] = pde[2] + recovery[2] + mtx[2] + js[2]

print("\n  L   J   finest    linf        rel linf    rms         "
      "total (h)  PDE (h)   recovery (h)  user (h)  sys (h)")
for L in LMAX:
    finest_mesh = '%d^2' % (base_mesh * 2 ** L) if base_mesh else '?'
    print("  %-3d %-3d %-9s %-11.4g %-11.4g %-11.4g %-10.3f %-9.3f %-13.3f %-9.3f %.3f"
          % (L, L - 2, finest_mesh, linferror[L], relinf[L], rmserror[L],
             t_total[L] / 3600, t_pde[L] / 3600, t_recovery[L] / 3600,
             t_user[L] / 3600, t_sys[L] / 3600))
print()

# --------------------------------------------- single level, if present ----
sl_linf, sl_total = {}, {}
for L in SL_LMAX:
    run = sl_results[L]
    err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
    sl_linf[L] = np.abs(err).max()
    pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
    sl_total[L] = pde[0] + recovery[0] + mtx[0] + js[0]

if SL_LMAX:
    print("  multilevel against single level, at matching target L")
    print("  L    ML linf      SL linf      ML time (h)  SL time (h)  speed-up")
    for L in SL_LMAX:
        speedup = sl_total[L] / t_total[L] if t_total[L] > 0 else float('nan')
        print("  %-4d %-12.4g %-12.4g %-12.3f %-12.3f %.1fx"
              % (L, linferror[L], sl_linf[L],
                 t_total[L] / 3600, sl_total[L] / 3600, speedup))
    print()

# ----------------------------------------------------------- figures -------
colours = plt.cm.rainbow(np.linspace(0, 1, max(len(LMAX), 4)))
markers = ["o", "v", "^", "<", ">", "s", "8"]
labels = ["L = %d (J = %d)" % (L, L - 2) for L in LMAX]


def _marker(i):
    return markers[i % len(markers)]


## 1. work-precision: error against the wall clock that bought it.
##
## With the single-level overlay this is the ML-vs-SL figure. Colour then
## encodes the METHOD and the point labels encode L -- one channel each. The
## per-L rainbow used otherwise would spend six colours on what the annotations
## already say, which is what made the overlay look unaffordable.
plt.figure()
if SL_LMAX:
    plt.plot([t_total[L] for L in LMAX], [linferror[L] for L in LMAX],
             color='tab:blue', marker='o', label='Multilevel')
    plt.plot([sl_total[L] for L in SL_LMAX], [sl_linf[L] for L in SL_LMAX],
             color='tab:red', marker='s', label='Single level')
    for L in LMAX:
        plt.annotate('L=%d' % L, (t_total[L], linferror[L]),
                     textcoords='offset points', xytext=(4, 4), fontsize=8)
    for L in SL_LMAX:
        plt.annotate('', xy=(t_total[L], linferror[L]),
                     xytext=(sl_total[L], sl_linf[L]),
                     arrowprops=dict(arrowstyle='->', color='0.6',
                                     linewidth=0.8, linestyle=':'))
    plt.legend()
else:
    handles = [plt.scatter(t_total[L], linferror[L], color=colours[i], marker=_marker(i))
               for i, L in enumerate(LMAX)]
    plt.legend(handles=handles, labels=labels)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wall clock computing time (s)')
plt.ylabel(r'$\ell_\infty$ error')
plt.tight_layout()
savefig("errorVsWorkWRTLs")

## 2. error against the finest grid the run actually reached.
if base_mesh:
    plt.figure()
    handles = [plt.scatter(base_mesh * 2 ** L, linferror[L], color=colours[i], marker=_marker(i))
               for i, L in enumerate(LMAX)]
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.xlabel('Finest grid of the run ($%d \\cdot 2^{L}$ points per side)' % base_mesh)
    plt.ylabel(r'$\ell_\infty$ error')
    plt.legend(handles=handles, labels=labels)
    plt.tight_layout()
    savefig("errorVsMeshWRTLs")

## 3. the headline: where the time goes, and the error, against L.
fig, ax1 = plt.subplots()
ax1.set_xlabel('L')
ax1.set_ylabel('time (h)')
for i, (series, label) in enumerate(((t_total, 'Total time'),
                                     (t_pde, 'PDE solve time'),
                                     (t_recovery, 'Sparse recovery time'))):
    ax1.plot(LMAX, [series[L] / 3600 for L in LMAX],
             color=colours[i], marker=_marker(i), label=label)
ax1.set_xticks(LMAX)

ax2 = ax1.twinx()
ax2.set_ylabel(r'$\ell_\infty$ error')
ax2.set_yscale('log')
ax2.plot(LMAX, [linferror[L] for L in LMAX],
         color=colours[3], marker=_marker(3), label='Final approximation error')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)
fig.tight_layout()
savefig("timesAndErrorWRTLs")

plt.show()
