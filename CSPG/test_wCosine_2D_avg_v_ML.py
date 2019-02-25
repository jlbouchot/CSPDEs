import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, WeightedCosine2D, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np
import argparse

# We still have to pass some inputs to the main:
# nb_max_iter: the number of max iteration for the greedy algo -> this is to be defined at the first place, it is not done at all for now. 
# recovery_algo: get to pick between the whtp (nb iter // norm residual // constance of the support), BPDN (norm on the residual), WOMP (norm residual? nb_iter? w(S^n) \leq sl?), wGHTP ('' '')
# ----> Make sure that all the arguments are passed accordingly. 


__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/02/24"


def get_sampling_type(sampling_name):
    switcher = {
        "pragmatic": WR.cs_pragmatic_m,
        "theoretic": WR.cs_theoretic_m,
        "new": WR.cs_theoretic_m_new,
		"p": WR.cs_pragmatic_m,
        "t": WR.cs_theoretic_m,
    }
    return switcher.get(sampling_name, WR.cs_pragmatic_m)



# def Main(outfile, d = 10, L_max = 4, orig_mesh_size = 2000):
def Main(outfile = "thatTest", d = 5, grid_points = tuple([2000,2000]), L_max = 4, algo_name = "whtp", gamma = 1.035, L_min = 1, sampling_name = "p", nb_iter = 500, epsilon = 1e-3, nb_tests = None, alpha = 2.0, abar = 4.3, imp = 1, w_cst = 0.5, dat_constant = 10, experiment_name = "weighted_cosine_avg_v_2D", tensor_based=True):


    		
    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModelML(WeightedCosine2D(d, alpha, imp, abar, w_cst), ConstantCoefficient(10.0),
                                       Average(), grid_points) 

	# Still have to concatenate the output file name with the parameters (i.e. d and h_0)
    test_result = outfile, None
    # test_result = '_'.join([algo_name, str(d), str(grid_points),outfile]), None
    for s in range(L_min,L_max+1,1): # s corresponds to the number of levels here
        ### Reconstruction Model
        v = np.hstack((np.repeat(gamma, d), [np.inf]))

        if tensor_based: 
            wr_model   = WR.WRModel(algo_name, WR.Operators.Cheb_Alt, v, 
                                get_sampling_type(sampling_name), WR.check_cs)
            prefix_npy = experiment_name + "Tensor_h"
        else: # The basic way.
            wr_model   = WR.WRModel(algo_name, WR.Operators.Chebyshev, v,
                                get_sampling_type(sampling_name), WR.check_cs)
            prefix_npy = experiment_name + "Classic_h"

		## Number of tests
        num_tests = nb_tests 

		### Execute test
        test_result = test(spde_model, wr_model, nb_iter, epsilon, s, [CrossCheck(num_tests)], dat_constant, prefix_npy + str(grid_points[0]) + "_", *test_result)
		## Don't forget to reset the original mesh
        spde_model.refine_mesh(2**(-(s-1)))


### Main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-d", "--nb-cosines", help="Number of random cosine and sine parameters", default=5, required=False)
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputDiffusionML", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-x", "--mesh-x", help="Size of the coarsest level (number of grid points) in the x direction", default=2000, required=False)
    parser.add_argument("-y", "--mesh-y", help="Size of the coarsest level (number of grid points) in the y direction", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=50, required=False)
    parser.add_argument("-e", "--tol-res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=1e-4, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant coefficients", default=1.035, required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default="pragmatic", required=False)
    parser.add_argument("-n", "--nb-tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-p", "--power", help="Power of the decay of the trigonometric expansion", default=2.0, required=False)
    parser.add_argument("-a", "--abar", help="Value of the mean field", default=4.3, required=False)
    parser.add_argument("-c", "--dat_constant", help="Multiplicative constant for the sparsity per level", default=10., required=False)
    parser.add_argument("-f", "--prefix-precompute", help="How should the precomputed data for this test be called?", default="", required=False)
    parser.add_argument("-b", "--better-compute", help="Should the computations be done on the fly, using tensor representation (Default is TRUE)", default="True", required=False)
    parser.add_argument("-i", "--fluctuation-importance", help="What is the importance of the fluctuations with respect to the mean field (default is 1)", default=1, required=False)
    parser.add_argument("-w", "--weight-cosine", help="How much weight the local cosine carries (Default = 0.5)", default=0.5, required=False)
    args = parser.parse_args()
	
    
    Main(args.output_file, int(args.nb_cosines), tuple([int(args.mesh_x),int(args.mesh_y)]), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), int(args.l_start), args.sampling, int(args.nb_iter), float(args.tol_res), None if args.nb_tests is None else int(args.nb_tests), float(args.power), float(args.abar), float(args.fluctuation_importance), float(args.weight_cosine), float(args.dat_constant), args.prefix_precompute, args.better_compute.lower()=="true")
    # Main(sys.argv[1])
