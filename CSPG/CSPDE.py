import numpy as np
import itertools
from collections import namedtuple


CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 'm', 'd', 'Z', 'y', 'A', 'w', 'result'])

def CSPDE(spde_model, wr_model, epsilon, s, cspde_result = None):
    ## 1. Create index set and draw random samples
    print("Generating J_s ...")

    # Compute "active index set" J_s
    J_s = J(s, wr_model.operator.theta, wr_model.weights)

    # Get total number of coefficients in tensorized chebyshev polynomial base
    N = len(J_s)

    # Calculate number of samples
    m = wr_model.get_m_from_s_N(s, N)

    # Get sample dimension
    d = len(J_s[0])

    if not cspde_result is None:
        assert d == cspde_result.d, "New sample space dimension is different from old sample space dimension."

    # Check whether this even an interesting case
    print("   It is N={0}, m={1} and d={2} ... ".format(N, m, d))
    wr_model.check(N, m)


    ## 2. Get samples from SPDEModel with given accuracy
    if cspde_result is None:
        # Draw m random samples in uniformly in [-1, 1] w.r.t. to the precondition measure associated to the given operator
        Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))

        print("\nComputing {0} SPDE sample approximations ...".format(m))
        # Get samples
        y = spde_model.samples(Z, epsilon)
    else:
        m_diff = m - cspde_result.m

        # Need to draw new samples ?
        if m_diff > 0:
            Z_new = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m_diff, d)))
            Z     = np.vstack((cspde_result.Z, Z_new))

            print("\nComputing {0} new SPDE sample approximations ...".format(m_diff))
            y     = np.hstack((cspde_result.y, spde_model.samples(Z_new, epsilon)))
        else:
            Z = cspde_result.Z[0:m]
            y = cspde_result.y[0:m]

    ## 3. Solve compressed sensing problem
    print("\nSolving compressed sensing problem ...")

    # Create sampling matrix and weights
    print("   Creating sample operator ...")
    A = wr_model.operator.create(J_s, Z)

    print("   Computing weights ...")
    w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), J_s)

    print("   Weighted minimization ...")
    result = wr_model.method(A, y, w, s, np.sqrt(m) * epsilon)


    return CSPDEResult(J_s, N, m, d, Z, y, 0, w, result) ## !!! REMOVED "A" FROM OUTPUT


def J(s, theta, v):
    print("s = {0}, theta = {1}, v = {2}".format(s,theta,v))
    # Function for generating all admissible indices over given index set S
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
    M = np.argmin(a <= A - T)
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


def calculate_weights(theta, v, J_s):
    return np.array([theta**np.count_nonzero(nu) * np.product(v[nu > 0]**nu[nu > 0]) for nu in J_s])
