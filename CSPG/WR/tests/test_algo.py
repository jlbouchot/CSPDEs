import pytest

import numpy as np
import WR
from .conftest import *



def tes_algorithms_Orig_unweighted(algo_orig, para_set_orig_unw):
    # x_truth = sparse_ground_truth(n)
    x_truth = para_set_orig_unw["x_truth"]
    # op = create_Operator(m,n)
    op = para_set_orig_unw["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set_orig_unw["y"]
    w = np.ones(para_set_orig_unw["n"])
    eta = para_set_orig_unw["epsilon"]
    s = para_set_orig_unw["s"]
    res = algo_orig(op, y, w, s, eta, para_set_orig_unw["maxiter"])
    # print("x_truth: ", x_truth)
    # print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10


def test_algorithms_Orig_weighted(algo_orig, para_set_orig_wei):
    # x_truth = sparse_ground_truth(n)
    x_truth = para_set_orig_wei["x_truth"]
    # op = create_Operator(m,n)
    op = para_set_orig_wei["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set_orig_wei["y"]
    # TODO: fix weight
    w = np.ones(para_set_orig_wei["n"])
    eta = para_set_orig_wei["epsilon"]
    s = para_set_orig_wei["s"]
    res = algo_orig(op, y, w, s, eta, para_set_orig_wei["maxiter"])
    # print("x_truth: ", x_truth)
    # print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-6


def test_algorithms_Alt_weighted(algo, para_set_Alt):
    # x_trut = sparse_ground_truth(n)
    x_truth = para_set_Alt["x_truth"]
    # op = create_Operator(m,n)
    op = para_set_Alt["operator"]
    # y = create_obs_y(op,x_truth)
    y = para_set_Alt["y"]
    w = para_set_Alt["weights"]
    # w = calculate_weights(para_set_Alt["theta"], para_set_Alt["v_weights"], para_set_Alt["J"])
    eta = para_set_Alt["epsilon"]
    # s = int(np.rint(para_set_Alt["density"]*para_set_Alt["n"]+0.01)) # 0.01 is an ofset because of rounding limitations in numpy.rint
    s = para_set_Alt["s"]
    res = algo(op, y, w, s, eta, para_set_Alt["maxiter"])
    # print("x_truth: ", x_truth)
    # print("res[\"x\"]", res[0])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-6
