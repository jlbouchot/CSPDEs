from num_int import vectorized_MC_integration as vmci

import numpy as np
import itertools

# vmci has the following signature: 
# vectorized_MC_integration(tot_tests, a_nu, batch_size = 100000, rhs = 1., variability = 1./6., mean_field = 5., regular_saves = True, dump_fname = None)
max_tests = int(100e6) # written this way for better readability
batch_size = int(1e6)
rhs = 1.0
variability = 1.0/6.0 # don't give too much importance to the local variation
mean_field = 5.0
dump_fname = "MCcoefficients_dim6_PWconstant_smallVar.npy"

max_degree = 2
d = 6

J = [list(x) for x in itertools.product(range(max_degree+1), repeat=d)]
for a_nu in J:
    vmci(max_tests, a_nu, batch_size, rhs, variability, mean_field, True, dump_fname)


