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
__lastmodified__ = "2026/10/28"

def wls_cgls(Operator, y, w, eta=1e-10, maxiter=200, print_every=False):
    ''' Matrix free weighted least squares solver using the conjugate gradient least squares (CGLS) method.
        Operator: linear operator with apply and apply_adj methods
        y: right-hand side vector
        w: weight vector (must be positive)
        eta: tolerance for the relative residual norm
        maxiter: maximum number of iterations
        print_every: print progress every 'print_every' iterations (if > 0)
    '''
    sw = np.sqrt(w)
    x = np.zeros(Operator.n)
    r = sw * y                                    # r = W^{1/2}(y - A x)
    s = Operator.apply_adj(sw * r)                # s = A* W (y - A x)
    p = s.copy()
    gamma = np.vdot(s, s).real
    gamma0 = gamma
    t = u.time_things()
    k = 0
    while k < maxiter and np.sqrt(gamma / gamma0) > eta: # and np.linalg.norm(r) > tol:
        q = sw * Operator.apply(p)
        alpha = gamma / np.vdot(q, q).real
        x += alpha * p
        r -= alpha * q
        s = Operator.apply_adj(sw * r)
        gamma_new = np.vdot(s, s).real
        p = s + (gamma_new / gamma) * p
        gamma = gamma_new
        k += 1
        if print_every and k % print_every == 0:
            print(f"WLS-CGLS iteration {k}: weighted residual = {np.linalg.norm(r):.3e}")
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
