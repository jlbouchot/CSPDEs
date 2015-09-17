import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, TrigCoefficient, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np


# We still have to pass some inputs to the main:
# d: dimension of the trigonometric expansion
# grid_points: the coarsest level of approximation
# nb_max_iter: the number of max iteration for the greedy algo -> this is to be defined at the first place, it is not done at all for now. 
# recovery_algo: get to pick between the whtp (nb iter // norm residual // constance of the support), BPDN (norm on the residual), WOMP (norm residual? nb_iter? w(S^n) \leq sl?), wGHTP ('' '')
# ----> Make sure that all the arguments are passed accordingly. 


def Main(outfile):
    ### Input parameters
    epsilon     = 1e-4 # Change from 8 to test on quinoa

    ## SPDEModel
    d     = 3 # Number of sines and cosines in the expansion
	
	# Corseast grid 
    grid_points = 200 # ~ 5e-3 

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModelML(TrigCoefficient(d, 1.0, 4.3), ConstantCoefficient(10.0),
                                       Average(), grid_points) 

	# Still have to concatenate the output file name with the parameters (i.e. d and h_0)
    test_result = '_'.join([str(d), str(grid_points),outfile]), None
    for s in range(1,3,1): # s corresponds to the number of levels here
        for gamma in np.arange(1.035, 1.04, 0.01)[::-1]:
            ### Reconstruction Model
            v = np.hstack((np.repeat(gamma, 2*d), [np.inf]))

            wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)

            ## Number of tests
            num_tests = 250 # change from 10 for Quinoa tests

            ### Execute test
            spde_model.refine_mesh(2**(-s))
            test_result = test(spde_model, wr_model, epsilon, s, [CrossCheck(num_tests)], *test_result)
            ## Don't forget to reset the orginal mesh


### Main
if __name__ == "__main__":
    Main(sys.argv[1])
