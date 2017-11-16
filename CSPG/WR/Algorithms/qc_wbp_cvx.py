import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

# cvx_path='~/install/lib/python2.7/site-packages/cvxopt-1.1.7-py2.7-linux-x86_64.egg'
# sys.path.append(cvx_path)
# import cvxopt as cvx
import cvxpy as cvx

from Result import *


def qc_wbp_cvx(Operator, y, w, s, eta, maxiter):
    x = cvx.Variable(Operator.n)

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) * x, 1)),
                [cvx.norm(Operator.A * x - y, 2) <= eta]).solve(solver=cvx.SCS) # Potentially try CVXOPT too

    return Result(np.array(x.value).flatten(), -1, 'Quadratically Constrained CVX')
