import pytest
import numpy as np
import WR
from scipy.sparse import rand
from WR.Operators.operator_from_matrix_Alt import univ_tensor_from_tensor_indices

parameter_set = {"epsilon": 1e-6,
                 "density": 0.05,
                 "n": 1000,
                 "maxiter": 100,
                 "d": 5,
                 "theta": np.sqrt(2),
                 "v" : 1.05
                 }

# >>> ws = np.array([1,1.4, 2.1,3,5])
# >>> ws
# array([1. , 1.4, 2.1, 3. , 5. ])
# >>> winv = ws**(-1)
# >>> winv
# array([1. , 0.71428571, 0.47619048, 0.33333333, 0.2 ])
# >>> cumDistribution = np.cumsum(winv)
# >>> cumDistribution
# array([1. , 1.71428571, 2.19047619, 2.52380952, 2.72380952])
# >>> cumDistribution = cumDistribution / cumDistribution[-1]
# >>> cumDistribution
# array([0.36713287, 0.62937063, 0.8041958 , 0.92657343, 1. ])
# >>> np.random.rand()
# 0.49981139949214004


# randstate = 42
# m = int(np.log(n) * density * n)
nb_test_runs = 6

def weighted_sparse_ground_truth(size, weights, density, randstate):
    np.random.seed(randstate)
    cumDist = np.cumsum(weights**(-1))
    cumDist = cumDist/cumDist[-1]
    truth = np.zeros(size)
    # sampling with replacement
    for i in range(0,int(size*density)):
        rand_loc = find_loc(np.random.rand(1)[0],cumDist) # because random.rand returns an array
        truth[rand_loc] = np.random.normal()
    return truth

def find_loc(x, dist):
    for i in range(0, len(dist)):
        if x < dist[i]:
            break
    return i

def sparse_ground_truth(size, density, randstate):
    matrix = rand(size, 1, density=density, format="csr", random_state=randstate).todense()
    return np.squeeze(np.asarray(matrix)) #squeeze needed because rand returns the deprecated matrix and we want to reduce dim

def create_Operator(m,n,randstate):
    np.random.seed(randstate)
    matrix = np.random.normal(size=(m,n))
    matrix = matrix/(np.sqrt(m))
    operator = WR.Operators.operator_from_matrix(WR.Operators.Chebyshev, matrix)
    return operator

def create_Operator_Alt(n,randstate):
    def base(x, k):
        return np.cos(k * np.arccos(x)) * np.sqrt(2)**(k>0)
    np.random.seed(randstate)
    J = calc_J(parameter_set["s"], parameter_set["theta"], parameter_set["v_weights"])
    m = WR.cs_theoretic_m_new(parameter_set["s"], len(J))
    d = parameter_set["d"]
    Z = WR.Operators.LD_bounded_operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
    operator = WR.Operators.operator_from_matrix_Alt(WR.Operators.Cheb_Alt, np.zeros((2,2)),
                                                     univ_tensor_from_tensor_indices(J, Z, base), J)
    return {"operator": operator, "J": J, "Z": Z, "m": m}

def create_obs_y(op,x):
    return op.apply(x)


@pytest.fixture(params=np.asarray(range(nb_test_runs))+14)
def para_set(request):
    parameter_set["randstate"] = request.param
    parameter_set["s"] = int(np.rint(parameter_set["density"]*parameter_set["n"]+0.01))
    # parameter_set["m"] = WR.cs_theoretic_m_new(parameter_set["s"], parameter_set["n"])
    parameter_set["m"] = int(np.ceil(2 * np.log10(parameter_set["s"])**2 * parameter_set["s"] * np.log10(parameter_set["n"])))
    # parameter_set["m"] = int(np.ceil(2 * np.log(parameter_set["s"])**2 * np.log(parameter_set["n"])))
    parameter_set["x_truth"] = sparse_ground_truth(parameter_set["n"],
                                                   parameter_set["density"],
                                                   parameter_set["randstate"])
    parameter_set["operator"] = create_Operator(parameter_set["m"],
                                                parameter_set["n"],
                                                parameter_set["randstate"])
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set


@pytest.fixture(params=np.asarray(range(nb_test_runs))+1)
def para_set_Alt(request):
    parameter_set["randstate"] = request.param
    parameter_set["v_weights"] = np.hstack((np.repeat(parameter_set["v"], parameter_set["d"]), [np.inf]))
    parameter_set["s"] = 10 #int(np.rint(parameter_set["density"]*parameter_set["n"]+0.01))
    # parameter_set["m"] = int(np.ceil(2 * np.log(parameter_set["s"])**2 * parameter_set["s"] * np.log(parameter_set["n"])))
    parameter_set["x_truth"] = sparse_ground_truth(parameter_set["n"],
                                                   parameter_set["density"],
                                                   parameter_set["randstate"])
    set_3 = create_Operator_Alt(parameter_set["n"],parameter_set["randstate"])
    parameter_set["J"] = set_3["J"]
    parameter_set["n"] = len(parameter_set["J"])
    parameter_set["m"] = set_3["m"]
    parameter_set["operator"] = set_3["operator"]
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set

@pytest.fixture(params=[WR.Algorithms.whtp, WR.Algorithms.wiht,
                        WR.Algorithms.womp])#, WR.Algorithms.exact_wbp_cvx,
                        # WR.Algorithms.qc_wbp_cvx])
#                         WR.Algorithms.Qc_wbp_precond_primaldual, WR.Algorithms.wcosamp,
#                         , , WR.Algorithms.primaldual])
def algo(request):
    return request.param

# needs to be imported from CSPDE_ML but we need to write an appropriate init file for CSPG to import from there
import itertools
def calc_J(s, theta, v):
    def iterate(M, a, B, S):
        def iterate_(B, S, p):
            L = []

            if len(S):
                while True:
                    # Take the highest index in S
                    r     = S[-1]

                    # Substract weight from B
                    B    -= a[r]

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
    A = np.log2(s/2.)
    a = 2 * np.log2(v)
    T = 2 * np.log2(theta)

    # Determine maximal M s.t. for j = 0 ... M-1 is a_j <= A - T
    # M is also the maximal support size
    # M = np.argmin(a <= A - T)
    M = np.argmin(a <= A )
    assert 0 != M, "Weight array too short. (Last element: {0}. Threshold: {1})".format(a[-1], A-T)

    # If A is non-negative the zero vector is always admissible
    assert A >= 0, "Negative A, i.e. sparsity less than 2."
    L = [np.zeros(M, dtype='int')]

    # Iterate through support sets of cardinality k = 1 ... M
    for k in range(1, M + 1):
        new_indices = []

        for S in itertools.combinations(range(M), k):
            new_indices += iterate(M, a, A - k*T, list(S))

        if [] == new_indices:
            break

        L += new_indices

    return L
