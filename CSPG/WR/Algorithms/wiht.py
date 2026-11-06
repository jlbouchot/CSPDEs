import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import *
from .weighted_quasi_abslargest import *
from .ls_solver import wls_cgls

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, RWTH Aachen, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2026/11/06"

def wiht(Operator, y, w, s, eta, maxiter, miniter, cs_ok=True, print_every=False, stag_tol=1e-3, stag_window=5):
    '''
    wiht implements a weighted iterative hard thresholding algorithm for sparse recovery. 
    In particular `wiht` solves a fixed point problem and is the `eta` is too restrictive, the algorithm may stagnate. 
    The parameters `stag_tol` and `stag_window` are used to detect stagnation and stop the algorithm early, preventing it from running to maxIter.
    '''
    if not cs_ok:
        ## TODO: ordinary least squares, not weighted. `wls_cgls` expects ROW
        ## weights (one per sample, length m), whereas `w` here is the
        ## length-N coefficient weight vector, which only steers the sparsity
        ## selection the N<m branch skips. No per-sample weights exist in this
        ## problem, so we pass ones. Revisit if heteroscedastic sample noise
        ## (e.g. level-dependent FE error) ever gives us genuine row weights.
        return wls_cgls(Operator, y, np.ones(np.size(y)), eta=1e-10, maxiter=maxiter, print_every=print_every)

    x = np.zeros(Operator.n)
    last_norm = 0
    k = 0

    t = u.time_things()
    residuum = y - Operator.apply(x)
    cur_norm = np.linalg.norm(Operator.apply(x) - y)

    while k < miniter or (k < maxiter and cur_norm > eta):
        x, dummy = weighted_quasi_abslargest_fast(x + Operator.apply_adj(residuum), 3 * s, w)
        k = k + 1
        residuum = y - Operator.apply(x)
        last_norm, cur_norm = cur_norm, np.linalg.norm(residuum)

        ## Stagnation guard
        if stag_window and np.isfinite(last_norm) and last_norm > 0 \
                and abs(last_norm - cur_norm) <= stag_tol * last_norm:
            n_stag += 1
        else:
            n_stag = 0
        if n_stag >= stag_window and k >= miniter:
            if print_every:
                print('WIHT stagnated after {0} iterations: {1} with less than {2:.3g} '
                      'relative change, residuum {3} against eta {4}.'
                      .format(k, stag_window, stag_tol, cur_norm, eta))
            break

        if print_every and print_every > 0:
            if k % print_every == 0:
                print('WIHT iteration {0}: norm of residuum = {1}'.format(k, cur_norm))
        if k > maxiter:
            print('WIHT did not converge after {0} iterations.'.format(k))
            break

    t = u.time_things(t)
    if print_every:
        print('WIHT finished after {0} iterations with norm of residuum = {1}'.format(k, cur_norm))

    # TODO: Better return code handling (tolerance, stagnation, maxiter reached, etc.)
    return Result(x, k, 'Weighted Iterative Hard Thresholding', cur_norm <= eta, t[0], t[1], t[2])
