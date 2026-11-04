import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import Result
from .weighted_quasi_abslargest import *
from .ls_solver import wls_cgls

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, RWTH Aachen, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2026/10/28"


def whtp(Operator, y, w, s, eta, maxiter, miniter, cs_ok=True, print_every=False):
    """Weighted hard thresholding pursuit.

    If the compressed-sensing condition is not satisfied, we fall back to a
    weighted least-squares CGLS solve. Otherwise we perform a thresholding-based
    support selection followed by a least-squares solve on that support.
    """
    if not cs_ok:
        return wls_cgls(Operator, y, w, eta=eta, maxiter=maxiter, print_every=print_every)

    x     = np.zeros(Operator.n)
    S_old = np.array([], dtype=int)
    k     = 0

    t = u.time_things()
    residuum = y - Operator.apply(x)
    cur_norm = np.linalg.norm(residuum)

    while k < miniter or (k < maxiter and cur_norm > eta):

        _, S = weighted_quasi_abslargest_fast(x + Operator.apply_adj(residuum), 3 * s, w)
        S = np.asarray(S, dtype=int)
        if S.size == S_old.size and np.array_equal(S, S_old):
            break

        ## All the following lines are equivalent to a weighted least-squares solve on the support S. Reuse the CGLS solver if you want to avoid the explicit matrix inversion.
        SubMatrix = Operator.genSubMatrix(S)
        x    = np.zeros(Operator.n)
        x[S] = np.linalg.lstsq(SubMatrix, y, rcond=None)[0]
        ## End of cgls.

        S_old = S
        k     = k + 1
        residuum = y - Operator.apply(x)
        cur_norm = np.linalg.norm(residuum)

        if print_every and print_every > 0:
            if k % print_every == 0:
                print('WHTP-new iteration {0}: norm of residual = {1}'.format(k, cur_norm))

    t = u.time_things(t)
    if print_every and not k >= maxiter:
        print('WHTP finished after {0} iterations with norm of residual = {1}'.format(k, cur_norm))
    if k >= maxiter:
        print('WHTP did not converge after {0} iterations.'.format(k))

    return Result(x, k, 'Weighted HTP (new)', k < maxiter, t[0], t[1], t[2])
