import numpy as np
import cvxpy as cvx

from Result import *


def qc_wbp_cvx(Operator, y, w, s, eta):
    x = cvx.Variable(Operator.n)

    cvx.Problem(cvx.Minimize(cvx.norm(np.diag(w) * x, 1)),
                [cvx.norm(Operator.A * x - y, 2) <= eta]).solve()

    return Result(np.array(x.value).flatten(), -1, 'Quadratically Constrained CVX')
