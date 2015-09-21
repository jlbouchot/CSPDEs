import numpy as np

from WR import eps
from Result import Result
from weighted_quasi_abslargest import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

def whtp(Operator, y, w, s, eta):
    x     = np.zeros(Operator.n)
    S_old = np.array([])
    k     = 0

    while True:
        if k > eta:
            print('WHTP did not converge after {0} iterations.'.format(k))
            break

        dummy, S = weighted_quasi_abslargest(x + Operator.apply_adj(y - Operator.apply(x)), 3 * s, w)

        if set(S) == set(S_old):
            break

        SubMatrix = Operator.genSubMatrix(S)

        x    = np.zeros(Operator.n)
        x[S] = np.linalg.lstsq(SubMatrix, y)[0]

        S_old = S
        k     = k + 1

    return Result(x, k, 'Weighted HTP')
