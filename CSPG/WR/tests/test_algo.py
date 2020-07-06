import pytest

import numpy as np
import WR
from WR.Algorithms.weighted_quasi_abslargest import *
from .conftest import *




def test_algorithms(algo, para_set):
    # x_truth = sparse_ground_truth(n)
    x_truth = para_set["x_truth"]
    # op = create_Operator(m,n)
    op = para_set["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set["y"]
    w = np.ones(parameter_set["n"])
    eta = para_set["epsilon"]
    s = int(np.rint(para_set["density"]*para_set["n"]+0.01)) # 0.01 is an ofset because of rounding limitations in numpy.rint
    res = WR.Algorithms.whtp(op, y, w, s, eta, para_set["maxiter"])
    print("x_truth: ", x_truth)
    print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10
