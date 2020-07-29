import pytest
import numpy as np
import WR
from scipy.sparse import rand
from WR.Operators.operator_from_matrix_Alt import univ_tensor_from_tensor_indices

parameter_set_test_config = {"epsilon": 1e-6,  # success in recovery ?
                             "density": 0.05,  # proportion of non zero components
                             "n": 1000,
                             # UNWEIGHTED ONLY Size of the ambient space, not in the parametric PDE like framework
                             "maxiter": 100,  # maximum number of iterations for the iterative algorithms
                             "d": 30,  # WEIGHTED ONLY number of _virtual_ operators
                             "theta": np.sqrt(2),
                             # WEIGHTED ONLY l-infinity norm of the univariate operators (see MLCSPG paper)
                             "v": 1.05,  # WEIGHTED ONLY Uniform weights for the operators
                             "s": 15  # will be OVERWRITTEN in UNWEIGHTED Case
                             }
nb_test_runs = 1  # number of tests


def weighted_sparse_ground_truth(size, weights, s, randstate):
    np.random.seed(randstate)  # set randomstate for reproducibility
    cumDist = np.cumsum(weights ** (-1))
    cumDist = cumDist / cumDist[-1]
    truth = np.zeros(size)  # initialize ground truth
    weighted_l0_norm = 0
    loc_set = set()  # set to remember the used locations
    while weighted_l0_norm < s:
        rand_loc = find_loc(np.random.rand(1)[0], cumDist)
        while rand_loc in loc_set:
            rand_loc = find_loc(np.random.rand(1)[0], cumDist)
        loc_set.add(rand_loc)
        curWeight = weights[rand_loc]
        print("Position {} with weight {}".format(rand_loc, curWeight))
        weighted_l0_norm = weighted_l0_norm + curWeight ** 2
        truth[rand_loc] = np.random.normal() / curWeight  # why divide through the weight?
        # truth[rand_loc] = 1.0/curWeight
    return truth


def find_loc(x, dist):
    for i in range(0, len(dist)):
        if x < dist[i]:
            break
    return i


# UNWEIGHTED ONLY
# def sparse_ground_truth(size, density, randstate):  # Only in the unweighted case
#     matrix = rand(size, 1, density=density, format="csr", random_state=randstate).todense()
#     return np.squeeze(
#         np.asarray(matrix))  # squeeze needed because rand returns the deprecated matrix and we want to reduce dim


# Original Operator
def create_Operator(m, n, randstate):  # Only in the unweighted case
    np.random.seed(randstate)
    # Create a random Gaussian matrix, which is known to fulfill (unweighted)RIP with high probability
    matrix = np.random.normal(size=(m, n))
    matrix = matrix / (np.sqrt(m))
    operator = WR.Operators.operator_from_matrix(WR.Operators.Chebyshev, matrix)
    return operator


# new Operator
def create_Operator_Alt(para_set):  # Only for the weighted case
    def base(x, k):
        return np.cos(k * np.arccos(x)) * np.sqrt(2) ** (k > 0)
    randstate = para_set["randstate"]
    np.random.seed(randstate)
    J = calc_J(para_set["s"], para_set["theta"], para_set["v_weights"])
    m = WR.cs_theoretic_m_new(para_set["s"], len(J))
    d = para_set["d"]
    Z = WR.Operators.LD_bounded_operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
    operator = WR.Operators.operator_from_matrix_Alt(WR.Operators.Cheb_Alt, np.zeros((2, 2)),
                                                     univ_tensor_from_tensor_indices(J, Z, base), J)
    return {"operator": operator, "J": J, "Z": Z, "m": m}


def create_obs_y(op, x):
    return op.apply(x)


@pytest.fixture(params=np.asarray(range(nb_test_runs)) + 42)
def para_set_orig_unw(request):
    parameter_set = dict(parameter_set_test_config)
    parameter_set["randstate"] = request.param
    parameter_set["v_weights"] = np.hstack((np.repeat(parameter_set["v"], parameter_set["d"]), [np.inf]))
    np.random.seed(parameter_set["randstate"])
    parameter_set["J"] = calc_J(parameter_set["s"], parameter_set["theta"], parameter_set["v_weights"])
    parameter_set["n"] = len(parameter_set["J"])
    parameter_set["m"] = WR.cs_theoretic_m_new(parameter_set["s"], len(parameter_set["J"]))
    print("n: {0}, m: {1}".format(parameter_set["n"], parameter_set["m"]))
    parameter_set["operator"] = create_Operator(parameter_set["m"],
                                                parameter_set["n"],
                                                parameter_set["randstate"])
    parameter_set["weights"] = np.ones(parameter_set["n"])
    parameter_set["x_truth"] = weighted_sparse_ground_truth(parameter_set["n"],
                                                            parameter_set["weights"],
                                                            parameter_set["s"],
                                                            parameter_set["randstate"])
    parameter_set["operator"] = create_Operator(parameter_set["m"],
                                                parameter_set["n"],
                                                parameter_set["randstate"])
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set

@pytest.fixture(params=np.asarray(range(nb_test_runs)) + 42)
def para_set_orig_wei(request):
    parameter_set = dict(parameter_set_test_config)
    parameter_set["randstate"] = request.param
    parameter_set["v_weights"] = np.hstack((np.repeat(parameter_set["v"], parameter_set["d"]), [np.inf]))
    np.random.seed(parameter_set["randstate"])
    parameter_set["J"] = calc_J(parameter_set["s"], parameter_set["theta"], parameter_set["v_weights"])
    parameter_set["n"] = len(parameter_set["J"])
    parameter_set["m"] = WR.cs_theoretic_m_new(parameter_set["s"], len(parameter_set["J"]))
    print("n: {0}, m: {1}".format(parameter_set["n"], parameter_set["m"]))
    parameter_set["operator"] = create_Operator(parameter_set["m"],
                                                parameter_set["n"],
                                                parameter_set["randstate"])
    parameter_set["weights"] = calculate_weights(parameter_set["theta"], parameter_set["v_weights"], parameter_set["J"])
    parameter_set["x_truth"] = weighted_sparse_ground_truth(parameter_set["n"],
                                                            parameter_set["weights"],
                                                            parameter_set["s"],
                                                            parameter_set["randstate"])
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set


@pytest.fixture(params=np.asarray(range(nb_test_runs)) + 42)
def para_set_Alt(request):
    parameter_set = dict(parameter_set_test_config)  # explicitly copy the basic test configuration from above
    # initialization of the parameters
    parameter_set["randstate"] = request.param
    parameter_set["v_weights"] = np.hstack((np.repeat(parameter_set["v"], parameter_set["d"]), [np.inf]))
    set_3 = create_Operator_Alt(parameter_set)
    parameter_set["J"] = set_3["J"]
    parameter_set["n"] = len(parameter_set["J"])
    parameter_set["m"] = set_3["m"]
    parameter_set["operator"] = set_3["operator"]
    parameter_set["weights"] = calculate_weights(parameter_set["theta"], parameter_set["v_weights"], parameter_set["J"])
    print("n: ", parameter_set["n"], " m: ", parameter_set["m"])
    parameter_set["x_truth"] = weighted_sparse_ground_truth(parameter_set["n"],
                                                            parameter_set["weights"],
                                                            parameter_set["s"],
                                                            parameter_set["randstate"])
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set


@pytest.fixture(params=[WR.Algorithms.whtp, WR.Algorithms.wiht,
                        WR.Algorithms.womp, WR.Algorithms.exact_wbp_cvx,
                        WR.Algorithms.qc_wbp_cvx, WR.Algorithms.wcosamp])
def algo_orig(request):
    return request.param


@pytest.fixture(params=[WR.Algorithms.whtp, WR.Algorithms.wiht,
                        WR.Algorithms.womp,
                        WR.Algorithms.wcosamp])
def algo(request):
    return request.param


# the following functions needs to be imported from CSPDE_ML but we need to write an appropriate init file for CSPG to import from there

def calculate_weights(theta, v, J_s):
    return np.array(
        [theta ** np.count_nonzero(nu) * np.product(v[np.where(nu > 0)] ** nu[np.where(nu > 0)]) for nu in J_s])


import itertools


def calc_J(s, theta, v):
    def iterate(M, a, B, S):
        def iterate_(B, S, p):
            L = []

            if len(S):
                while True:
                    # Take the highest index in S
                    r = S[-1]

                    # Substract weight from B
                    B -= a[r]

                    # Increase nu_r by one
                    p[r] += 1

                    if len(S) > 1:
                        # If there is more than one index left, recurse with remaining indices and 'remaining weight' B
                        e = iterate_(B, S[0:-1], p.copy())
                    else:
                        # If B - sum over a_j is larger than 0, add this multiindex
                        if B >= 0:
                            e = [p.copy()]
                        else:
                            e = []

                    # If the list of new multiindices is empty, there is nothing left to be done
                    # thanks to the monotonicity of a
                    if not e:
                        break

                    # Add found multiindices
                    L += e

            return L

        return iterate_(B, S, np.zeros(M, dtype='int'))

    # Set A and a as in Theorem 5.2
    A = np.log2(s / 2.)
    a = 2 * np.log2(v)
    T = 2 * np.log2(theta)

    # Determine maximal M s.t. for j = 0 ... M-1 is a_j <= A - T
    # M is also the maximal support size
    # M = np.argmin(a <= A - T)
    M = np.argmin(a <= A)
    assert 0 != M, "Weight array too short. (Last element: {0}. Threshold: {1})".format(a[-1], A - T)

    # If A is non-negative the zero vector is always admissible
    assert A >= 0, "Negative A, i.e. sparsity less than 2."
    L = [np.zeros(M, dtype='int')]

    # Iterate through support sets of cardinality k = 1 ... M
    for k in range(1, M + 1):
        new_indices = []

        for S in itertools.combinations(range(M), k):
            new_indices += iterate(M, a, A - k * T, list(S))

        if [] == new_indices:
            break

        L += new_indices

    return L
