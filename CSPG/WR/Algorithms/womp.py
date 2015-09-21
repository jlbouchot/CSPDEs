import numpy as np

from WR import eps
from Result import Result

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

def womp(Operator, y, w, s, eta):
    x = np.zeros(Operator.n)
    S = []
    k = 0
    s = 3*s # This should allow for ample error (13s being the lowest upper bound)

    while s > 0 and np.linalg.norm(Operator.apply(x) - y) > eta:
        j    = np.argmax(np.abs(Operator.apply_adj(y - Operator.apply(x))) * 1./w)

        if j in S:
            break

        S.append(j)

        x[S] = np.linalg.lstsq(Operator.genSubMatrix(S), y)[0]

        s    = s - 1
        k    = k + 1

    return Result(x, k, 'Weighted OMP')
