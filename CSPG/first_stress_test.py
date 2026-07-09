# Test file to make sure things work with the MLCSPG and fenics/dolfin configuration
import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, WeightedCosine2D, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np
import argparse


__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/07/06"


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
def Main():
# def Main(outfile = "thatTest", epsilon = 1e-3, experiment_name = "weighted_cosine_avg_p_2D"):

    experiment_name = "stresstest"
    nb_tests = 20
    dict_config = {'d': 5, 'J': 2, "L": 4, "h0": tuple([20, 20]), "vj": 1.035, "weightCosine": 0.5, "nbSamples": "t", "Tensor": True, 't': 1, "tprime": 1, 'p0': 0.3, "p": 0.3, "s_J": 5, "s_L": 10, "trig_power": 2., "abar": 4.3, "energy_fluctuations": 1., "algo": "wiht", "iter": 20, "tolres": 1e-3, "ansatz": 0, "no_compute": True, "alpha": 0.25}

    # Adapt to the first approximating level (via a single level approach)
    grid_points = tuple(int(2**(2)*dummy) for dummy in dict_config["h0"])

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModelML(WeightedCosine2D(dict_config["d"], dict_config["trig_power"], dict_config["energy_fluctuations"], dict_config["abar"], dict_config["weightCosine"]), ConstantCoefficient(10.0),
                                       Average(), grid_points) 

    # Still have to concatenate the output file name with the parameters (i.e. d and h_0)
    outfile = "noComputeTestWHTP"
    test_result = None
    # test_result = outfile, None
    v = np.hstack((dict_config["vj"]*np.power([val+1 for val in range(dict_config["d"]) for dummy_variable in (0,1)], dict_config["alpha"]), [np.inf]))

    wr_model   = WR.WRModel("wiht", WR.Operators.Cheb_Alt, v, 
                            get_sampling_type(dict_config["nbSamples"]), WR.check_cs)

	## Number of tests
    num_tests = nb_tests 

	### Execute test
#     test_result = test(spde_model, wr_model, nb_iter, epsilon, L_min, L_max, [CrossCheck(num_tests)], dat_constant, p, p0, t_0, t_prime, const_sJ, ansatz_space, prefix_npy + str(grid_points[0]) + "_", *test_result)
    test_result = test(spde_model, wr_model, dict_config, [CrossCheck(num_tests)], experiment_name, test_result)


    dict_config["algo"] = "bpdn"
    wr_model   = WR.WRModel("bpdn", WR.Operators.Chebyshev, v,
                            get_sampling_type(dict_config["nbSamples"]), WR.check_cs)
    outfile = "doComputeTestCVX"
    test_result = None
    test_result = test(spde_model, wr_model, dict_config, [CrossCheck(num_tests)], experiment_name, *test_result)


### Main
if __name__ == "__main__":

    # Main(args.output_file, int(args.nb_cosines), tuple([int(args.mesh_x),int(args.mesh_y)]), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), int(args.l_start), args.sampling, int(args.nb_iter), float(args.tol_res), None if args.nb_tests is None else int(args.nb_tests), float(args.power), float(args.abar), float(args.fluctuation_importance), float(args.weight_cosine), float(args.dat_constant), args.prefix_precompute, args.better_compute.lower()=="true", int(args.ansatz_space), float(args.t_0), float(args.t_prime), float(args.smooth_0), float(args.smooth_t), float(args.const_sJ), False if args.nb_tests is None else args.no_compute, float(args.exponent))
    # Main(sys.argv[1])
    Main()
