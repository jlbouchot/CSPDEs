import numpy as np
import itertools
from collections import namedtuple

import time


__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery'])

def CSPDE_ML(spde_model, wr_model, unscaledNbIter, epsilon, L=1, dat_constant = 10, cspde_result = None):
	lvl_by_lvl_result = []
	for oneLvl in xrange(0,L):
		# sl = 10+np.max([2**(L-oneLvl),1])

		sl = np.floor(dat_constant*2**(L-oneLvl))
		print("Computing level {0} from a total of {1}. Current sparsity = {2}".format(oneLvl+1,L,sl))
		## 1. Create index set and draw random samples
		print("Generating J_s ...")
	
		# Compute "active index set" J_s
		J_s = J(sl, wr_model.operator.theta, wr_model.weights)

		# Get total number of coefficients in tensorized chebyshev polynomial base
		N = len(J_s)
	
		# Calculate number of samples
		m = wr_model.get_m_from_s_N(sl, N)
	
		# Get sample dimension
		d = len(J_s[0])
	
		# if not cspde_result is None:
		#     assert d == cspde_result.d, "New sample space dimension is different from old sample space dimension."
	
		# Check whether this even an interesting case
		print("   It is N={0}, m={1} and d={2} ... ".format(N, m, d))
		wr_model.check(N, m)
	
		Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
		print("\nComputing {0} SPDE sample approximations ...".format(m))
		# Get samples
		t_start = time.time()
		if oneLvl != 0:
			y_old = spde_model.samples(Z)
			spde_model.refine_mesh()
		else:
			y_old = np.zeros(m)
		y_new = spde_model.samples(Z)
		t_stop = time.time()
		t_samples = t_stop-t_start

		## 3. Solve compressed sensing problem
		print("\nSolving compressed sensing problem ...")
	
		# Create sampling matrix and weights
		print("   Creating sample operator ...")
		t_start = time.time()
		A = wr_model.operator.create(J_s, Z)
                t_stop = time.time()
		t_matrix = t_stop-t_start
		
		print("   Computing weights ...")
		w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), J_s)
	
		print("   Weighted minimization ...")
		t_start = time.time()
		result = wr_model.method(A, y_new-y_old, w, sl, epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
		t_stop = time.time()
		t_recovery = t_stop-t_start
		# result = wr_model.method(A, y_new-y_old, w, sl, np.sqrt(m) *epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
		lvl_by_lvl_result.append(CSPDEResult(J_s, N, sl, m, d, Z, y_new-y_old, 0, w, result, t_samples, t_matrix, t_recovery))
		print("\n\tRecovery time: {0} \t Building the Matrix: {1} \t Computing the samples {2}\n".format(t_recovery, t_matrix, t_samples))
	
	
	return lvl_by_lvl_result


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


def calculate_weights(theta, v, J_s):
    return np.array([theta**np.count_nonzero(nu) * np.product(v[np.where(nu > 0)]**nu[np.where(nu > 0)]) for nu in J_s])
    # return np.array([theta**np.count_nonzero(nu) * np.product(v[nu > 0]**nu[nu > 0]) for nu in J_s]) old version for earlier python distributions
