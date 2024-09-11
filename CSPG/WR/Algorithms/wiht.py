import numpy as np

from WR import eps
from .Result import *
from .weighted_quasi_abslargest import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2024, INRIA, RWTH Aachen, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2024/09/11"

def wiht(Operator, y, w, s, eta, maxiter):
    x         = np.zeros(Operator.n)
    last_norm = 0
    k         = 0

    while np.linalg.norm(Operator.apply(x) - y) > eta:
        residuum = y - Operator.apply(x)
        cur_norm = np.linalg.norm(residuum)

        x, dummy  = weighted_quasi_abslargest(x + Operator.apply_adj(residuum), 3 * s, w)
        last_norm = cur_norm
        k         = k + 1

        if k > maxiter:
            print('WIHT did not converge after {0} iterations.'.format(k))
            break
    print("Weighted Recovery:WIHT stopped after {} iterations. Current residual is {}".format(k, np.linalg.norm(Operator.apply(x) - y)))

    return Result(x, k, 'Weighted Iterative Hard Thresholding')
