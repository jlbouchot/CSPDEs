import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from WR import eps
from .Result import *
from .weighted_quasi_abslargest import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, RWTH Aachen, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2026/08/03"

def wiht(Operator, y, w, s, eta, maxiter, print_every=False):
    x         = np.zeros(Operator.n)
    last_norm = 0
    k         = 0

    t = u.time_things()
    residuum = y - Operator.apply(x)
    cur_norm = np.linalg.norm(Operator.apply(x) - y)

    while cur_norm > eta:
        
        x, dummy  = weighted_quasi_abslargest(x + Operator.apply_adj(residuum), 3 * s, w)
        last_norm = cur_norm
        k         = k + 1
        residuum = y - Operator.apply(x)
        cur_norm = np.linalg.norm(residuum)

        if print_every and print_every > 0:
            if k % print_every == 0:
                print('WIHT iteration {0}: norm of residuum = {1}'.format(k, cur_norm))
        if k > maxiter:
            print('WIHT did not converge after {0} iterations.'.format(k))
            break
    
    t = u.time_things(t)
    if print_every:
        print('WIHT finished after {0} iterations with norm of residuum = {1}'.format(k, cur_norm))

    return Result(x, k, 'Weighted Iterative Hard Thresholding', k <= maxiter, t[0], t[1], t[2])
