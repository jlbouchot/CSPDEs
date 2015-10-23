import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML, LinearCoefficient, ConstantCoefficient, Integration, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np

import math

def Main(outfile = "testPiecewiseConstantDiffusion", grid_points = 2000, L_max = 3, algo_name = "whtp", gamma = 1.035, abar = 5, variability = None, L_min = 1):
    
	
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
		
    ## SPDEModel
    d  = 13
    epsilon = 50 # Is used for the number of iterations in whtp

    mesh_size = int(d*math.floor(float(grid_points)/float(d)))

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    if variability is None:
        variability = abar/(d+1)
    #upper_bound = abar + variability
    #a  = [LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
    #      LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
    #      LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound)] # Make sure you have d of those
    spde_model = PiecewiseConstantDiffusionFEMModelML(abar, ConstantCoefficient(1.0), variability, mesh_size, Average())

    test_result = outfile, None

    for s in range(L_min,L_max+1,1):
        ## Reconstruction Model
        v = [gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, gamma, np.inf] # This has to be done better too

        wr_model = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                              WR.cs_pragmatic_m, WR.check_cs)

        num_tests = 2000 # change from 10 for Quinoa tests

		## Don't forget to reset the original mesh
        spde_model.refine_mesh(2**(-s))

        ### Execute test
        test_result = test(spde_model, wr_model, epsilon, s, [CrossCheck(num_tests)], *test_result)


### Main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputPiecewiseConstantDiffusionML", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-m", "--mesh-size", help="Size of the coarsest level (number of grid points)", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=500, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant weights", default=1.035, required=False)
    parser.add_argument("-a", "--abar", help="Mean (constant) diffusion field", default=5, required=False)
    parser.add_argument("-v", "--variability", help="Variations (uncertainty) in the mean field per section - None is equivalent to abar/(nb sections + 1)", default=None, required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)

    args = parser.parse_args()
    
    Main(args.output_file, int(args.mesh_size), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), float(args.abar), None if args.variability is None else float(args.variability), int(args.l_start))
