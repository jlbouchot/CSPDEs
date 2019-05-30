import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, CosineCoef3D, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np
import argparse

# We still have to pass some inputs to the main:
# nb_max_iter: the number of max iteration for the greedy algo -> this is to be defined at the first place, it is not done at all for now. 
# recovery_algo: get to pick between the whtp (nb iter // norm residual // constance of the support), BPDN (norm on the residual), WOMP (norm residual? nb_iter? w(S^n) \leq sl?), wGHTP ('' '')
# ----> Make sure that all the arguments are passed accordingly. 


__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/05/30"


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
def Main(outfile = "testCosine3D", d = 5, grid_points = tuple([20, 20, 20]), L_max = 3, algo_name = "whtp", gamma = 1.035, L_min = 3, sampling_name = "p", nb_iter = 50, epsilon = 1e-3, nb_tests = None, alpha = 2.0, abar = 4.3, imp = 1,  dat_constant = 10, experiment_name = "avg_v_3D", tensor_based=True, ansatz_space = 0, t_0 = 1, t_prime = 1, p0 = 1./4., p = 3./10., const_sJ = 5):

    
    # Adapt to the first approximating level (via a single level approach)
    grid_points = tuple(int(2**(L_min-1)*dummy) for dummy in grid_points)
		
    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModelML(CosineCoef3D(d, alpha, imp, abar), ConstantCoefficient(10.0),
                                       Average(), grid_points) 

	# Still have to concatenate the output file name with the parameters (i.e. d and h_0)
    test_result = outfile, None
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
        num_tests = nb_tests # change from 10 for Quinoa tests


		### Execute test
        test_result = test(spde_model, wr_model, nb_iter, epsilon, L_min, L_max, [CrossCheck(num_tests)], dat_constant, p, p0, t_0, t_prime, const_sJ, ansatz_space, prefix_npy + str(grid_points[0]) + "_", *test_result)
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
    parser.add_argument("-z", "--mesh-z", help="Size of the coarsest level (number of grid points) in the z direction", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=50, required=False)
    parser.add_argument("-e", "--tol-res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=1e-4, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant coefficients", default=1.035, required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default="pragmatic", required=False)
    parser.add_argument("-n", "--nb-tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-p", "--power", help="Power of the decay of the trigonometric expansion", default=2.0, required=False)
    parser.add_argument("-a", "--abar", help="Value of the mean field", default=4.3, required=False)
    parser.add_argument("-c", "--dat-constant", help="Multiplicative constant for the sparsity per level", default=10., required=False)
    parser.add_argument("-f", "--prefix-precompute", help="How should the precomputed data for this test be called?", default="", required=False)
    parser.add_argument("-b", "--better-compute", help="Should the computations be done on the fly, using tensor representation (Default is TRUE)", default="True", required=False)
    parser.add_argument("-i", "--fluctuation-importance", help="What is the importance of the fluctuations with respect to the mean field (default is 1)", default=1, required=False)
    parser.add_argument("-j", "--ansatz-space", help="What type of Ansatz space is used? (Default is 0)", default="0", required=False)
    parser.add_argument("--t_0", help="What is the smoothness of the data (Default is 1)", default="1", required=False)
    parser.add_argument("--t_prime", help="What is the smoothness of the functional (Default is 1)", default="1", required=False)
    parser.add_argument("--smooth_0", help="What kind of smoothness in the original space can be expected (Default is 1/4)", default="0.25", required=False)
    parser.add_argument("--smooth_t", help="What kind of smoothness in the smooth space can be expected (Default is 3/10)", default="0.3", required=False)
    parser.add_argument("--const_sJ", help="What is the expected constant in the expression of s_J (Default is 10)", default="10", required=False)
    args = parser.parse_args()
	
    
    Main(args.output_file, int(args.nb_cosines), tuple([int(args.mesh_x),int(args.mesh_y), int(args.mesh_z)]), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), int(args.l_start), args.sampling, int(args.nb_iter), float(args.tol_res), None if args.nb_tests is None else int(args.nb_tests), float(args.power), float(args.abar), float(args.fluctuation_importance), float(args.dat_constant), args.prefix_precompute, args.better_compute.lower()=="true",int(args.ansatz_space), float(args.t_0), float(args.t_prime), float(args.smooth_0), float(args.smooth_t), float(args.const_sJ))
    # Main(sys.argv[1])
