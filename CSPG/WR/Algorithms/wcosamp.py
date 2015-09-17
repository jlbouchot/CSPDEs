import numpy as np

from WR import eps
from Result import Result
from weighted_quasi_abslargest import *


def wcosamp(Operator, y, w, s, eta):
    x = np.zeros(Operator.n)
    u = np.zeros(Operator.n)
    U_old = []

    last_norm = np.inf
    k         = 0

    while k < s:
        if k > 100:
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

    return Result(x, k, 'Weighted CoSaMP')
