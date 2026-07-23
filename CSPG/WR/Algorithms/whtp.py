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

def whtp(Operator, y, w, s, eta, maxiter):
    x     = np.zeros(Operator.n)
    S_old = np.array([])
    k     = 0

    t = u.time_things()

    while True:
        if k > maxiter:
            print('WHTP did not converge after {0} iterations.'.format(k))
            break

        dummy, S = weighted_quasi_abslargest(x + Operator.apply_adj(y - Operator.apply(x)), 3 * s, w)
        if set(S) == set(S_old) or np.linalg.norm(Operator.apply(x) - y)/np.linalg.norm(y) <= eta:
            print('WHTP Converged after {0} iterations. Norm of residual {1}'.format(k,np.linalg.norm(Operator.apply(x) - y)))
            break

        SubMatrix = Operator.genSubMatrix(S)

        x    = np.zeros(Operator.n)
        x[S] = np.linalg.lstsq(SubMatrix, y)[0]

        S_old = S
        k     = k + 1

    t = u.time_things(t)

    return Result(x, k, 'Weighted HTP', k<=maxiter, t[0], t[1], t[2])
