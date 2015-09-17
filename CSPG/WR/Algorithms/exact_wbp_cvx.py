import numpy as np
import cvxpy as cvx

from WR import eps
from Result import *


def exact_wbp_cvx(Operator, y, w, s, eta):
    x = cvx.Variable(Operator.n)

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) * x, 1)),
                [cvx.norm(Operator.A * x - y, 1) <= eps]).solve()

    return Result(np.array(x.value).flatten(), -1, 'Exact CVX')
