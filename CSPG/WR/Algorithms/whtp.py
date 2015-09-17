import numpy as np

from WR import eps
from Result import Result
from weighted_quasi_abslargest import *


def whtp(Operator, y, w, s, eta):
    x     = np.zeros(Operator.n)
    S_old = np.array([])
    k     = 0

    while True:
        if k > 200:
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
