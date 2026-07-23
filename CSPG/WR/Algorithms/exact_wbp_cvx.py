import numpy as np
## Add utilities to the path

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/07/21"

import cvxpy as cvx

from WR import eps
from .Result import *


def exact_wbp_cvx(Operator, y, w, s, eta, maxiter):
    x = cvx.Variable(Operator.n)

    t = u.time_things()

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) @ x, 1)),
                [cvx.norm(Operator.A @ x - y, 1) <= eps]).solve(solver=cvx.SCS) # Potentially try CVXOPT too

    t = u.time_things(t)

    return Result(np.array(x.value).flatten(), -1, 'Exact CVX', True, t[0], t[1], t[2])
