from compare_true_coefs import get_computed_coefs as gcc
from num_int import coefs_MC_wrapper as cmw
import numpy as np
import math 

infile = 'smallDim6_PWLD_whtp_var_small_grid2000_L3_g105_datconstant20'
outfile = 'testingMCint_smallvar_datconstant20'

estimated_coefs, J = gcc(infile,outfile, J = None)

nb_tests = [10000, 100000, 1000000, 10000000]

print '\t\t Comparison MC vs l1 recovery\n'
print '\t nu \t\tMC 10000 \t\tMC 100000\t\tMC 1000000\t\tMC 10000000\t\tl1value \t\t Abs error'

for idx_nu, one_nu in enumerate(J): 
    if np.sum(one_nu) < 4: # We'll display only the nus with relatively low degrees
        values = np.zeros(len(nb_tests))
        for idx, nmc in enumerate(nb_tests):
            values[idx] = cmw(nmc, one_nu, False, 1.0, 1./6., 5.)
        print '{}\t\t {:.5f}\t\t {:.5f}\t\t {:.5f}\t\t {:.5f}\t\t {:.5f}\t\t {}'.format(one_nu,values[0], values[1], values[2], values[3], estimated_coefs[idx_nu], abs(values[3]-estimated_coefs[idx_nu]))
