"""Exp2: approximation error and compute time as the dimension d grows.

Reads one experiment directory of ``d_<nb_cosines>`` shelves and plots them
against the ground truth. Unlike the other figures this needs ONE reference per
dimension -- d is part of the identity, so every point has its own folder.

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
from utilities.utils import getComputeTimes, ground_truth_for, load_runs, resolve_qoi

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

# The sweep, as nb_cosines. nb dimensions = 2*nb_cosines.
NB_COSINES = [3, 5, 7, 10, 13, 16, 20, 24, 28]      # testing/launchExp2.sh's grid
RESULTS = os.path.join(TESTING, 'results', 'avgqoi_dgrows_preQoIfix')

NB_TESTS = 2000

# Configuration for the ground truth.
GT_CFG = os.path.join(HERE, '..', 'data', 'Exp2_dvaries.ini')
GT_SEED = 20261106
#Mesh type used for computing the ground truth.
GT_MESH = 'rich640-1280-2560'
# ============================================================================


# ------------------------------------------------------------------ load ---
runs = load_runs(RESULTS, dict((k, 'd_%d' % k) for k in NB_COSINES))
# Check whether the ground truth and runs agree on the QOI.
# One could check more things to make sure everything agrees...
resolve_qoi(runs)

# ------------------------------------------- errors and timings, per d -----
linferror, rmserror, relinf = {}, {}, {}
t_total, t_pde, t_recovery = {}, {}, {}
t_user, t_sys = {}, {}

for k in NB_COSINES:
    run = runs[k]
    d = 2 * k
    ## One reference per dimension: d is part of the identity hash, so each of
    ## these resolves to its own GTresults_d<d>H<mesh>_<hash> folder.
    y_GT, Z, identity = ground_truth_for({k: run}, GT_CFG, GT_MESH, GT_SEED, d,
                                         NB_TESTS, base_dir=TESTING)
    err = np.asarray(run.wr_model.estimate_ML_samples(run.cspde_result, Z)).ravel() - y_GT
    linferror[k] = np.abs(err).max()
    rmserror[k] = np.linalg.norm(err) / np.sqrt(err.size)
    relinf[k] = linferror[k] / np.abs(y_GT).max()

    pde, recovery, mtx, js = getComputeTimes(run.cspde_result)
    ## Each is [wall clock, user, system]; seconds throughout.
    t_total[k] = pde[0] + recovery[0] + mtx[0] + js[0]
    t_pde[k] = pde[0]
    t_recovery[k] = recovery[0]
    t_user[k] = pde[1] + recovery[1] + mtx[1] + js[1]
    t_sys[k] = pde[2] + recovery[2] + mtx[2] + js[2]

print("\n  nb_cos  d    linf        rel linf    rms         "
      "total (h)  PDE (h)   recovery (h)  user (h)  sys (h)")
for k in NB_COSINES:
    print("  %-7d %-4d %-11.4g %-11.4g %-11.4g %-10.3f %-9.3f %-13.3f %-9.3f %.3f"
          % (k, 2 * k, linferror[k], relinf[k], rmserror[k],
             t_total[k] / 3600, t_pde[k] / 3600, t_recovery[k] / 3600,
             t_user[k] / 3600, t_sys[k] / 3600))
print()

# ----------------------------------------------------------- figures -------
ds = [2 * k for k in NB_COSINES]
colours = plt.cm.rainbow(np.linspace(0, 1, 4))

# Verification of the growth of the compute time with respect to the number of dimensions.
fig, ax1 = plt.subplots()
ax1.set_xlabel('d')
ax1.set_ylabel('time (h)')
for i, (series, label) in enumerate(((t_total, 'Total time'),
                                     (t_pde, 'PDE solve time'),
                                     (t_recovery, 'Sparse recovery time'))):
    ax1.plot(ds, [series[k] / 3600 for k in NB_COSINES],
             color=colours[i], marker='o', label=label)
ax1.set_xticks(ds)

ax2 = ax1.twinx()
ax2.set_ylabel(r'$\ell_\infty$ error')
ax2.set_yscale('log')
ax2.plot(ds, [linferror[k] for k in NB_COSINES],
         color=colours[3], marker='<', label='Final approximation error')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)
fig.tight_layout()
plt.savefig("timesAndErrorWRTDs.eps")
plt.savefig("timesAndErrorWRTDs.jpg")

# Where the time goes, as a share of each run -- the shape of the growth.
plt.figure()
bottom = np.zeros(len(ds))
for i, (series, label) in enumerate(((t_pde, 'PDE solves'),
                                     (t_recovery, 'Sparse recovery'))):
    share = np.array([series[k] / t_total[k] for k in NB_COSINES])
    plt.bar(ds, share, bottom=bottom, width=1.5, color=colours[i], label=label)
    bottom += share
plt.bar(ds, 1.0 - bottom, bottom=bottom, width=1.5, color=colours[2],
        label='Matrix + index set')
plt.xlabel('d')
plt.ylabel('share of the total wall clock')
plt.xticks(ds)
plt.ylim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig("timeSharesWRTDs.eps")
plt.savefig("timeSharesWRTDs.jpg")

plt.show()
