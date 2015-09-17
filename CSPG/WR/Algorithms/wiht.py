import numpy as np

from WR import eps
from Result import *
from weighted_quasi_abslargest import *


def wiht(Operator, y, w, s, eta):
    x         = np.zeros(Operator.n)
    last_norm = 0
    k         = 0

    while np.linalg.norm(Operator.apply(x) - y) > eta:
        residuum = y - Operator.apply(x)
        cur_norm = np.linalg.norm(residuum)

        x, dummy  = weighted_quasi_abslargest(x + Operator.apply_adj(residuum), 3 * s, w)
        last_norm = cur_norm
        k         = k + 1

        if k > 1000:
            print('WIHT did not converge after {0} iterations.'.format(k))
            break

    return Result(x, k, 'Weighted Iterative Hard Thresholding')
