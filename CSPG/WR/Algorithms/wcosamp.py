import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import Result
from .weighted_quasi_abslargest import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2026/07/21"

def wcosamp(Operator, y, w, s, eta, maxiter):
    x = np.zeros(Operator.n)
    u = np.zeros(Operator.n)
    U_old = []

    last_norm = np.inf
    k         = 0

    t = u.time_things()

    while k < s:
        if k > maxiter:
            print('WCoSaMP did not converge after {0} iterations.'.format(k))
            break

        x, S = weighted_quasi_abslargest(u, 3*s, w)
        U = list(set(S).union(weighted_quasi_abslargest(Operator.apply_adj(y - Operator.apply(x)), 2*(3*s), w)[1]))

        if set(U) == set(U_old):
            break

        SubMatrix = Operator.genSubMatrix(U)

        u    = np.zeros(Operator.n)
        u[U] = np.linalg.lstsq(SubMatrix, y)[0]

        U_old = U
        k     = k + 1
    t = u.time_things(t)

    return Result(x, k, 'Weighted CoSaMP', k <= maxiter, t[0], t[1], t[2])
