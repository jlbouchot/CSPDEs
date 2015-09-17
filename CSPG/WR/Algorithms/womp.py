import numpy as np

from WR import eps
from Result import Result


def womp(Operator, y, w, s, eta):
    x = np.zeros(Operator.n)
    S = []
    k = 0

    while s > 0 and np.linalg.norm(Operator.apply(x) - y) > eta:
        j    = np.argmax(np.abs(Operator.apply_adj(y - Operator.apply(x))) * 1./w)

        if j in S:
            break

        S.append(j)

        x[S] = np.linalg.lstsq(Operator.genSubMatrix(S), y)[0]

        s    = s - 1
        k    = k + 1

    return Result(x, k, 'Weighted OMP')
