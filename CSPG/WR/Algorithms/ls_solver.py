import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import *
from .weighted_quasi_abslargest import *

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, RWTH Aachen, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/28"
__lastmodified__ = "2026/11/10"

def wls_cgls(Operator, y, w, eta=1e-10, maxiter=200, print_every=False, tol_res=None):
    ''' Matrix free weighted least squares solver using the conjugate gradient least squares (CGLS) method.
        Operator: linear operator with apply and apply_adj methods
        y: right-hand side vector
        w: ROW weight vector, one entry per sample (must be positive)
        eta: RELATIVE tolerance on the normal-equation residual
             ||A* W (y - A x)||, measured against its value at x = 0. It is a
             convergence test for the least-squares problem itself and carries
             no units of the data.
        tol_res: ABSOLUTE tolerance on the weighted residual
             ||W^(1/2) (y - A x)||_2. This is the SAME quantity the thresholding
             algorithms stop on, on the same scale as the per-level
             epsilon = tol_res * sqrt(m) that CSPDE_v2.run_plan derives, so a
             level that has spent its error budget stops here exactly as it
             would under wiht/whtp. None (the default) disables it and leaves a
             pure least-squares solve. Whichever of the two criteria is met
             first ends the iteration.
        maxiter: maximum number of iterations
        print_every: print progress every 'print_every' iterations (if > 0)

        NOTE: the residual tested against tol_res is the recursively updated one
        (r -= alpha q), not a freshly evaluated y - A x. Over the few tens of
        iterations this solver needs, the drift is orders of magnitude below any
        tol_res the schedule produces; recompute it explicitly if maxiter is
        ever raised into the thousands.
    '''
    sw = np.sqrt(w)
    x = np.zeros(Operator.n)
    r = sw * y                                    # r = W^{1/2}(y - A x)
    s = Operator.apply_adj(sw * r)                # s = A* W (y - A x)
    p = s.copy()
    gamma = np.vdot(s, s).real
    gamma0 = gamma
    rnorm = np.linalg.norm(r)                     # ||W^(1/2)(y - A x)||_2 at x = 0
    t = u.time_things()
    k = 0
    while k < maxiter and np.sqrt(gamma / gamma0) > eta \
            and (tol_res is None or rnorm > tol_res):
        q = sw * Operator.apply(p)
        alpha = gamma / np.vdot(q, q).real
        x += alpha * p
        r -= alpha * q
        rnorm = np.linalg.norm(r)
        s = Operator.apply_adj(sw * r)
        gamma_new = np.vdot(s, s).real
        p = s + (gamma_new / gamma) * p
        gamma = gamma_new
        k += 1
        if print_every and k % print_every == 0:
            print(f"WLS-CGLS iteration {k}: weighted residual = {rnorm:.3e}")
    if print_every:
        why = ("residual budget tol_res = %.3e reached" % tol_res) \
            if (tol_res is not None and rnorm <= tol_res) \
            else ("maxiter" if k >= maxiter else "least-squares convergence")
        print(f"WLS-CGLS finished after {k} iterations, weighted residual = "
              f"{rnorm:.3e} ({why})")
    t = u.time_things(t)
    return Result(x, k, 'Weighted least squares (CGLS)', k < maxiter, t[0], t[1], t[2])

def wls_direct(A, y, w, method="qr"):
    ''' Direct solver for weighted least squares problems using either Cholesky or QR decomposition.
        A: matrix of the linear system
        y: right-hand side vector
        w: weight vector (must be positive)
        method: "chol" for Cholesky, "qr" for QR decomposition
    '''
    sw = np.sqrt(w)
    Aw = A * sw[:, None]          # O(mn), no explicit diag(W) matrix
    yw = y * sw
    if method == "chol":          # cheapest, squares the condition number
        G = Aw.conj().T @ Aw      # n x n Gram matrix
        return cho_solve(cho_factor(G), Aw.conj().T @ yw)
    Q, R = qr(Aw, mode="economic")  # stable default
    return solve_triangular(R, Q.conj().T @ yw)
