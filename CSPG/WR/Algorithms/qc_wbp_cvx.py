import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

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

import cvxpy as cvx

from .Result import *


def qc_wbp_cvx(Operator, y, w, s, eta, maxiter):
    t = u.time_things()
    x = cvx.Variable(Operator.n)

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) @ x, 1)),
                [cvx.norm(Operator.A @ x - y, 2) <= eta]).solve(solver=cvx.SCS) # Potentially try CVXOPT too
    t = u.time_things(t)

    return Result(np.array(x.value).flatten(), -1, 'Quadratically Constrained CVX', True, t[0], t[1], t[2])
