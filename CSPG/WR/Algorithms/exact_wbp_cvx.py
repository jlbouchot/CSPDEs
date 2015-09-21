import numpy as np
import sys

# cvx_path='~/install/lib/python2.7/site-packages/cvxopt-1.1.7-py2.7-linux-x86_64.egg'
# sys.path.append(cvx_path)
# import cvxopt as cvx
import cvxpy as cvx

from WR import eps
from Result import *


def exact_wbp_cvx(Operator, y, w, s, eta):
    x = cvx.Variable(Operator.n)

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) * x, 1)),
                [cvx.norm(Operator.A * x - y, 1) <= eps]).solve()

    return Result(np.array(x.value).flatten(), -1, 'Exact CVX')
