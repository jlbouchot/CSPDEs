from compare_true_coefs import get_computed_coefs as gcc
from num_int import coefs_MC_wrapper as cmw
import numpy as np

infile = 'smallDim6_PWLD_whtp_var_big_grid2000_L4_g105'
outfile = 'testingMCint'

estimated_coefs, J = gcc(infile,outfile, J = None)

nb_tests = [10000, 100000, 500000, 1000000]

print '\t\t Comparison MC vs l1 recovery\n'
print '\t nu \tMC 10000 \tMC 100000\t500000\t1000000\tl1value\n'

for idx_nu, one_nu in enumerate(J): 
    if np.max(one_nu) < 2: # We'll display only the nus with relatively low degrees
        values = zeors(len(nb_tests))
        for idx, nmc in enumerate(nb_tests):
            values[idx] = cmw(nmc, one_nu)
        print '{} {:.4f}  {:.4f}  {:.4f}  {:.4f}  {:.4f}'.format(one_nu,values[0], values[1], values[2], values[3], estimated_coefs[idx_nu])
