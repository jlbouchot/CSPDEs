import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, TrigCoefficient, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np
import argparse

# We still have to pass some inputs to the main:
# nb_max_iter: the number of max iteration for the greedy algo -> this is to be defined at the first place, it is not done at all for now. 
# recovery_algo: get to pick between the whtp (nb iter // norm residual // constance of the support), BPDN (norm on the residual), WOMP (norm residual? nb_iter? w(S^n) \leq sl?), wGHTP ('' '')
# ----> Make sure that all the arguments are passed accordingly. 


__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"


# def Main(outfile, d = 10, L_max = 4, orig_mesh_size = 2000):
def Main():
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-d", "--nb-cosines", help="Number of random cosine and sine parameters", default=5, required=False)
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputDiffusionML", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-m", "--mesh-size", help="Size of the coarsest level (number of grid points)", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=500, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)

    args = parser.parse_args()
	
	# Corseast grid 
    outfile = args.output_file
    grid_points = int(args.mesh_size)
    d = int(args.nb_cosines)
    L_max = int(args.nb_level)
    algo_name = args.recovery_algo.lower()
    if algo_name == 'whtp': # Really have to find a way to deal with the epsilon/eta/nbIter parameter
        epsilon = 50 # This will be rescaled later
    elif algo_name == 'wiht':
	    epsilon = 1e-4 
    elif algo_name == 'womp':
	    epsilon = 1e-4 
    elif algo_name == 'bpdn':
	    epsilon = 1e-4 
    else: 
        epsilon = 50 # This will be rescaled later
		
    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModelML(TrigCoefficient(d, 1.0, 4.3), ConstantCoefficient(10.0),
                                       Average(), grid_points) 

	# Still have to concatenate the output file name with the parameters (i.e. d and h_0)
    test_result = '_'.join([algo_name, str(d), str(grid_points),outfile]), None
    for s in range(1,L_max+1,1): # s corresponds to the number of levels here
        for gamma in np.arange(1.055, 1.06, 0.01)[::-1]:
            ### Reconstruction Model
            v = np.hstack((np.repeat(gamma, 2*d), [np.inf]))

            wr_model   = WR.WRModel(algo_name, WR.Operators.Chebyshev, v,
                                    WR.cs_theoretic_m_new, WR.check_cs)
            #wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
            #                        WR.cs_pragmatic_m, WR.check_cs) # or cs_theoretic_m

			## Number of tests
            num_tests = 1000 # change from 10 for Quinoa tests

			## Don't forget to reset the original mesh
            spde_model.refine_mesh(2**(-s))
			### Execute test
            test_result = test(spde_model, wr_model, epsilon, s, [CrossCheck(num_tests)], *test_result)


### Main
if __name__ == "__main__":
    Main()
    # Main(sys.argv[1])