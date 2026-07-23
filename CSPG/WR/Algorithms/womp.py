import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import Result

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

def womp(Operator, y, w, s, eta, maxiter):
    x = np.zeros(Operator.n)
    S = []
    k = 0
    s = 13*s # This is the result that ensure the best 's' sparse recovery according to T. Zhang

    t = u.time_things()
    while s > 0 and np.linalg.norm(Operator.apply(x) - y) > eta:
        print('Norm of the error {0}, while the objective norm is {1} and the ratio is {2}'.format(np.linalg.norm(Operator.apply(x) - y),np.linalg.norm(y),np.linalg.norm(Operator.apply(x) - y)/np.linalg.norm(y)))
        j    = np.argmax(np.abs(Operator.apply_adj(y - Operator.apply(x))) * 1./w)

        if j in S:
            break

        S.append(j)

        x[S] = np.linalg.lstsq(Operator.genSubMatrix(S), y)[0]

        s    = s - 1
        k    = k + 1

    t = u.time_things(t)
    return Result(x, k, 'Weighted OMP', k <= maxiter, t[0], t[1], t[2])
