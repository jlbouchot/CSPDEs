import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML, LinearCoefficient, ConstantCoefficient, Integration, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np


def Main(outfile):

    ## SPDEModel
    d  = 9
    epsilon = 50 # Is used for the number of iterations in whtp

    mesh_size = d*200 # corresponds to 1800 grid points for our example

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    abar = 5.0
    variability = 1.0/(d+1)
    upper_bound = abar + variability
    a  = [LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
          LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
          LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound)] # Make sure you have d of those
    spde_model = PiecewiseConstantDiffusionFEMModelML(a, ConstantCoefficient(1.0), mesh_size, Average())

    test_result = outfile, None

    for s in range(1,3,1):
        for gamma in np.arange(1.035, 1.05, 0.01)[::-1]:
            ## Reconstruction Model
            v = [gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, np.inf] # This has to be done better too

            wr_model = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                  WR.cs_pragmatic_m, WR.check_cs)

								  
								  
            num_tests = 1000 # change from 10 for Quinoa tests

			## Don't forget to reset the original mesh
            spde_model.refine_mesh(2**(-s))
            
            ### Execute test
            test_result = test(spde_model, wr_model, epsilon, s, [CrossCheck(num_tests)], *test_result)


### Main
if __name__ == "__main__":
    Main(sys.argv[1])
