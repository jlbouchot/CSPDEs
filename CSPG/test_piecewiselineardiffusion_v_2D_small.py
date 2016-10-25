import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML2D_small, LinearCoefficient, ConstantCoefficient, Integration, Average

import argparse

from Check_ML import test, CrossCheck

import sys
import numpy as np

import math

def get_sampling_type(sampling_name):
    switcher = {
        "pragmatic": WR.cs_pragmatic_m,
        "theoretic": WR.cs_theoretic_m,
        "new": WR.cs_theoretic_m_new,
		"p": WR.cs_pragmatic_m,
        "t": WR.cs_theoretic_m,
    }
    return switcher.get(sampling_name, WR.cs_pragmatic_m)

def Main(outfile = "testPiecewiseConstantDiffusion2D", grid_points = tuple([200,200]), L_max = 3, algo_name = "whtp", gamma = 1.035, abar = 5, variability = None, L_min = 1, sampling_name = "p", nb_iter = 50, epsilon = 1e-3, nb_tests = None, dat_constant = 10):
    
	
    # if algo_name == 'whtp': # Really have to find a way to deal with the epsilon/eta/nbIter parameter
        # epsilon = 1e-4 # This will be rescaled later
    # elif algo_name == 'wiht':
	    # epsilon = 1e-4 
    # elif algo_name == 'womp':
	    # epsilon = 1e-4 
    # elif algo_name == 'bpdn':
	    # epsilon = 1e-4 
    # else: 
        # epsilon = 1e-4 # This will be rescaled later
	
	
    ## SPDEModel
    d  = 16

    # mesh_size = int(d*math.floor(float(grid_points)/float(d)))

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    if variability is None:
        variability = abar/(d+1)
    #upper_bound = abar + variability
    #a  = [LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
    #      LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
    #      LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound)] # Make sure you have d of those
    spde_model = PiecewiseConstantDiffusionFEMModelML2D_small(abar, ConstantCoefficient(1.0), variability, grid_points, Average())

    test_result = outfile, None

    for s in range(L_min,L_max+1,1):
        ## Reconstruction Model
        v = np.hstack((np.repeat(gamma, d), [np.inf])) # The 2D PW constant has been implemented as a 5 times 5 grid

        wr_model = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v, get_sampling_type(sampling_name), WR.check_cs)

        num_tests = nb_tests # change from 10 for Quinoa tests

		## Don't forget to reset the original mesh
        spde_model.refine_mesh(2**(-s))

        ### Execute test
        test_result = test(spde_model, wr_model, nb_iter, epsilon, s, [CrossCheck(num_tests)], dat_constant, *test_result)


### Main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputPiecewiseConstantDiffusionML", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-x", "--mesh-x", help="Size of the coarsest level (number of grid points) in the x direction", default=2000, required=False)
    parser.add_argument("-y", "--mesh-y", help="Size of the coarsest level (number of grid points) in the y direction", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=500, required=False)
    parser.add_argument("-e", "--tol-res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=1e-4, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant weights", default=1.035, required=False)
    parser.add_argument("-a", "--abar", help="Mean (constant) diffusion field", default=5, required=False)
    parser.add_argument("-v", "--variability", help="Variations (uncertainty) in the mean field per section - None is equivalent to abar/(nb sections + 1)", default=None, required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)
    parser.add_argument("-n", "--nb-tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default="pragmatic", required=False)
    parser.add_argument("-c", "--dat-constant", help="Multiplicative constant for the sparsity per level", default=10., required=False)

    args = parser.parse_args()
    
    Main(args.output_file, tuple([int(args.mesh_x),int(args.mesh_x)]), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), float(args.abar), None if args.variability is None else float(args.variability), int(args.l_start), args.sampling, int(args.nb_iter), float(args.tol_res), None if args.nb_tests is None else float(args.nb_tests), args.dat_constant)

