import pytest

import numpy as np
import WR
from .conftest import *




def tes_algorithms(algo, para_set):
    # x_truth = sparse_ground_truth(n)
    x_truth = para_set["x_truth"]
    # op = create_Operator(m,n)
    op = para_set["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set["y"]
    w = np.ones(parameter_set["n"])
    eta = para_set["epsilon"]
    s = para_set["s"]
    res = algo(op, y, w, s, eta, para_set["maxiter"])
    # print("x_truth: ", x_truth)
    # print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10

def test_algorithms_Alt(algo, para_set_Alt):
    # x_truth = sparse_ground_truth(n)
    x_truth = para_set_Alt["x_truth"]
    # op = create_Operator(m,n)
    op = para_set_Alt["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set_Alt["y"]
    w = calculate_weights(para_set_Alt["theta"], para_set_Alt["v_weights"], para_set_Alt["J"])
    print("w")
    # print(5/0)
    eta = para_set_Alt["epsilon"]
    s = int(np.rint(para_set_Alt["density"]*para_set_Alt["n"]+0.01)) # 0.01 is an ofset because of rounding limitations in numpy.rint
    res = WR.Algorithms.whtp(op, y, w, s, eta, para_set_Alt["maxiter"])
    print("x_truth: ", x_truth)
    print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10

# needs to be imported from CSPDE_ML for testing reasons copied here
def calculate_weights(theta, v, J_s):
    return np.array([theta**np.count_nonzero(nu) * np.product(v[np.where(nu > 0)]**nu[np.where(nu > 0)]) for nu in J_s])
