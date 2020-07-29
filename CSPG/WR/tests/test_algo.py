import numpy as np


def test_algorithms_Orig_unweighted(algo_orig, para_set_orig_unw):
    x_truth = para_set_orig_unw["x_truth"]
    op = para_set_orig_unw["operator"]
    y = para_set_orig_unw["y"]
    w = np.ones(para_set_orig_unw["n"])
    w = para_set_orig_unw["weights"]
    eta = para_set_orig_unw["epsilon"]
    s = para_set_orig_unw["s"]
    res = algo_orig(op, y, w, s, eta, para_set_orig_unw["maxiter"])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10


def test_algorithms_Orig_weighted(algo_orig, para_set_orig_wei):
    x_truth = para_set_orig_wei["x_truth"]
    op = para_set_orig_wei["operator"]
    y = para_set_orig_wei["y"]
    w = para_set_orig_wei["weights"]
    eta = para_set_orig_wei["epsilon"]
    s = para_set_orig_wei["s"]
    res = algo_orig(op, y, w, s, eta, para_set_orig_wei["maxiter"])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10


def test_algorithms_Alt_weighted(algo, para_set_Alt):
    x_truth = para_set_Alt["x_truth"]
    op = para_set_Alt["operator"]
    y = para_set_Alt["y"]
    w = para_set_Alt["weights"]
    eta = para_set_Alt["epsilon"]
    s = para_set_Alt["s"]
    res = algo(op, y, w, s, eta, para_set_Alt["maxiter"])
    assert np.linalg.norm(res[0]-x_truth) <= 1e-10
