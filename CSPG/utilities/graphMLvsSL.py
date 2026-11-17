"""Multilevel against single level, at matching target levels L.

Reads two experiment directories -- the multilevel runs (``L_val_<L>``) and the
single-level runs (``SL_val_<L>``, from testing/launchSL.sh) -- which share one
configuration and therefore one reference. The point of the figure is the
work-precision plane: for the method to pay for itself the multilevel points
must sit below and to the left of the single-level ones.

Both sets are checked against the reference's QoI together, so a directory pair
that disagrees is refused rather than plotted: scoring against the wrong
functional yields a smooth, plausible curve instead of an error.
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
__lastmodified__ = "2026/11/15"


# ============================== settings ====================================
HERE = os.path.dirname(os.path.abspath(__file__))
TESTING = os.path.join(HERE, '..', 'testing')

LMAX = [2, 3, 4, 5, 6]
ML_RESULTS = os.path.join(TESTING, 'results', 'avgqoi_h0_nbDisc_Fixed_preQoIfix')
SL_RESULTS = os.path.join(TESTING, 'results', 'avgqoi_SingleLevel')

NB_TESTS = 1000

# GT configurations
GT_CFG = os.path.join(HERE, '..', 'data', 'Exp1_hfvaries.ini')
GT_SEED = 20261106
GT_MESH = 'rich640-1280-2560'
# ============================================================================


if not os.path.isdir(SL_RESULTS):
    raise SystemExit(
        "no single-level results in %s.\n"
        "Generate them first:  cd testing && bash launchSL.sh\n"
        "(same .ini as the multilevel runs, so they share one reference)."
        % SL_RESULTS)

# ------------------------------------------------------------------ load ---
ml = load_runs(ML_RESULTS, dict((L, 'L_val_%d' % L) for L in LMAX))
sl = load_runs(SL_RESULTS, dict((L, 'SL_val_%d' % L) for L in LMAX))

d = ml[LMAX[-1]].cspde_result[0].d

# -------------------------------------------------------- ground truth -----
both = {}
for L in LMAX:
    both[('ML', L)] = ml[L]
    both[('SL', L)] = sl[L]
y_GT, Z, identity = ground_truth_for(both, GT_CFG, GT_MESH, GT_SEED, d,
                                     NB_TESTS, base_dir=TESTING)

# ------------------------------------------------- errors and timings ------
def score(runs):
    linf, total = {}, {}
    for L in LMAX:
        run = runs[L]
        err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
        linf[L] = np.abs(err).max()
        pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
        total[L] = pde[0] + recovery[0] + mtx[0] + js[0]   # [0] is wall clock
    return linf, total


ml_linf, ml_time = score(ml)
sl_linf, sl_time = score(sl)

print("\n  L    ML linf      SL linf      ML time (h)  SL time (h)  speed-up")
for L in LMAX:
    speedup = sl_time[L] / ml_time[L] if ml_time[L] > 0 else float('nan')
    print("  %-4d %-12.4g %-12.4g %-12.3f %-12.3f %.1fx"
          % (L, ml_linf[L], sl_linf[L], ml_time[L] / 3600, sl_time[L] / 3600, speedup))
print()

# ----------------------------------------------------------- figures -------
## The work-precision plane: the statement the comparison is actually making.
plt.figure()
plt.plot([ml_time[L] for L in LMAX], [ml_linf[L] for L in LMAX],
         color='tab:blue', marker='o', label='Multilevel')
plt.plot([sl_time[L] for L in LMAX], [sl_linf[L] for L in LMAX],
         color='tab:red', marker='s', label='Single level')
for L in LMAX:
    plt.annotate('L=%d' % L, (ml_time[L], ml_linf[L]),
                 textcoords='offset points', xytext=(4, 4), fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wall clock computing time (s)')
plt.ylabel(r'$\ell_\infty$ error')
plt.legend()
plt.tight_layout()
savefig("errorVsWorkMLvsSL")

## Cost at each target level, which is where the speed-up is read off.
plt.figure()
plt.plot(LMAX, [ml_time[L] / 3600 for L in LMAX],
         color='tab:blue', marker='o', label='Multilevel')
plt.plot(LMAX, [sl_time[L] / 3600 for L in LMAX],
         color='tab:red', marker='s', label='Single level')
plt.yscale('log')
plt.xticks(LMAX)
plt.xlabel('L')
plt.ylabel('Wall clock time (h)')
plt.legend()
plt.tight_layout()
savefig("timesMLvsSL")

plt.show()
