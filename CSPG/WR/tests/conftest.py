import pytest
import numpy as np
import WR
from scipy.sparse import rand
from WR.Operators.operator_from_matrix_Alt import univ_tensor_from_tensor_indices

parameter_set = {"epsilon": 1e-6, # success in recovery ?
                 "density": 0.05, # proportion of non zero components
                 "n": 1000, # UNWEIGHTED ONLY Size of the ambient space, not in the parametric PDE like framework
                 "maxiter": 100, # maximum number of iterations for the iterative algorithms
                 "d": 30, # WEIGHTED ONLY number of _virtual_ operators
                 "theta": np.sqrt(2), # WEIGHTED ONLY l-infinity norm of the univariate operators (see MLCSPG paper)
                 "v" : 1.05 # WEIGHTED ONLY Uniform weights for the operators
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

def weighted_sparse_ground_truth(size, weights, s, randstate):
    np.random.seed(randstate)
    cumDist = np.cumsum(weights**(-1))
    cumDist = cumDist/cumDist[-1]
    truth = np.zeros(size)
    # sampling with replacement
    weighted_l0_norm = 0
    while weighted_l0_norm < s: 
        rand_loc = find_loc(np.random.rand(1)[0],cumDist)
        curWeight = weights[rand_loc]
        print("Position {} with weight {}".format(rand_loc, curWeight))
        weighted_l0_norm = weighted_l0_norm + curWeight**2
        truth[rand_loc] = np.random.normal()/curWeight
        truth[rand_loc] = 1.0/curWeight
    # for i in range(0,int(size*density)):
    #     rand_loc = find_loc(np.random.rand(1)[0],cumDist) # because random.rand returns an array
    #     truth[rand_loc] = np.random.normal()/weights[rand_loc]
    #     truth[rand_loc] = 1.0/weights[rand_loc]
    return truth

def find_loc(x, dist):
    for i in range(0, len(dist)):
        if x < dist[i]:
            break
    return i

def calculate_weights(theta, v, J_s):
    return np.array([theta**np.count_nonzero(nu) * np.product(v[np.where(nu > 0)]**nu[np.where(nu > 0)]) for nu in J_s])

#only useful in the unweighted case
def sparse_ground_truth(size, density, randstate): # Only in the unweighted case
    matrix = rand(size, 1, density=density, format="csr", random_state=randstate).todense()
    return np.squeeze(np.asarray(matrix)) #squeeze needed because rand returns the deprecated matrix and we want to reduce dim

def create_Operator(m,n,randstate): # Only in the unweighted case
    np.random.seed(randstate)
    # Create a random Gaussian matrix, which is known to fulfill (unweighted)RIP with high probability
    matrix = np.random.normal(size=(m,n))
    matrix = matrix/(np.sqrt(m))
    operator = WR.Operators.operator_from_matrix(WR.Operators.Chebyshev, matrix)
    return operator

def create_Operator_Alt(randstate): # Only for the weighted case 
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
    # initialization of the parameters including the ones already fixed at the top of the file
    parameter_set["randstate"] = request.param
    parameter_set["v_weights"] = np.hstack((np.repeat(parameter_set["v"], parameter_set["d"]), [np.inf]))
    # fixed the sparsity parameter s at the moment
    # while s=10 the resulting number of len(J) is 56
    # I don't understand how the parameter s correlates to the density value
    parameter_set["s"] = 20 #int(np.rint(parameter_set["density"]*parameter_set["n"]+0.01))
    # parameter_set["m"] = int(np.ceil(2 * np.log(parameter_set["s"])**2 * parameter_set["s"] * np.log(parameter_set["n"])))
    # to create the operator we don't need the parameter n (which is the size of the ground truth) why don't I need it
    # or where do I have to insert it in the function create_Operator_Alt()?
    set_3 = create_Operator_Alt(parameter_set["randstate"])
    parameter_set["J"] = set_3["J"]
    # the problem why the first test did return a segmentation fault and the following test did not was:
    # here I reset the parameter n based on the size of J
    # for the first test then x_truth was initialized for the size 1000
    # then after initialization of x we set it to len(J)
    # which lead it to initialize x_truth in the secound run with len(J) instead of 1000
    parameter_set["n"] = len(parameter_set["J"])
    # with this configuration all the test seem to work yet it still holds that m>>n with 427 > 56
    # I don't know how to properly set or redefine this WR.cs_theoretic_m_new for a good result (m<n)
    parameter_set["m"] = set_3["m"]
    parameter_set["operator"] = set_3["operator"]
    parameter_set["weights"] = calculate_weights(parameter_set["theta"], parameter_set["v_weights"], parameter_set["J"])
    parameter_set["true_dim"] = 1000
    # now I moved the calculation of the ground truth after the initialization process of J so that I have the correct value for n
    # and instead of the dimension 1000 we are now working on the dimension 56 I think which is why I am confused
    print("n: ", parameter_set["n"], " m: ", parameter_set["m"])
    parameter_set["x_truth"] = weighted_sparse_ground_truth(parameter_set["n"],
                                                    parameter_set["weights"],
                                                    parameter_set["s"],
                                                   parameter_set["randstate"])
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
