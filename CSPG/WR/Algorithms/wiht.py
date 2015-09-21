import numpy as np

from WR import eps
from Result import *
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
