import MLCSPG_param_parser as MLCSPGParse

import WR

import utils as u

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, WeightedCosine2D, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np
import argparse


__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019-2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/01/22"


def Main(main_cfg, pde_cfg, sparse_cfg):
    '''
    Main function for calling 
    '''

    print(40 * "=")
    print((5 * " ") + " MAIN CONFIGURATION")
    for (k,v) in main_cfg.items():
        print((8 * " ") + k + " : " + str(v))
    print(40 * "=")
    print((10 * " ") + " PDE SOLVER CONFIGURATION")
    for (k,v) in pde_cfg.items():
        print((8 * " ") + k + " : " + str(v))
    print(40 * "=")
    print((5 * " ") + " SPARSE SOLVER CONFIGURATION")
    for (k,v) in sparse_cfg.items():
        print((8 * " ") + k + " : " + str(v))

    print("\n" + "\t" + 3*"-" + " LAUNCHING BASIC TESTCASE " + 3*"-")

    # Adapt to the first approximating level (via a single level approach)
    grid_points = tuple(int(2**(main_cfg["l_start"])*dummy) for dummy in [main_cfg["mesh_x"], main_cfg["mesh_y"], main_cfg["mesh_z"]][0:main_cfg["n"]])

    # TODO: Handle parametrization of the linear coefficients in the operator expansion better
    # Suggestion: Consider a default case of WeightedCosine with a set of default coefs
    # Then, potentially, add yet another config for this with the appropriate parameters
    coef_expansion = WeightedCosine2D(main_cfg["nb_cosines"], main_cfg["power"], main_cfg["fluctuation_importance"], main_cfg["abar"], main_cfg["weight_cosine"]) # THIS IS FOR 2D things! NEeds to be updated with some utils
    spde_model = DiffusionFEMModelML(coef_expansion, ConstantCoefficient(10.0),
                                       Average(), grid_points) 
    
    test_result = None
        ### Reconstruction Model
    # WHY DID WE HAVE THIS LINE LIKE THIS ??? v = np.hstack((gamma*np.power([val+1 for val in range(d) for dummy_variable in (0,1)], exponent), [np.inf]))
    v = np.hstack((main_cfg["gamma"]*np.power([val+1 for val in range(main_cfg["nb_cosines"])], main_cfg["exponent"]), [np.inf]))

    if main_cfg["do_tensor"]: 
        wr_model   = WR.WRModel(sparse_cfg["recovery_algo"], WR.Operators.Cheb_Alt, v, 
                            u.get_sampling_type(main_cfg["sampling"]), WR.check_cs)
    else: # The basic way.
        wr_model   = WR.WRModel(sparse_cfg["recovery_algo"], WR.Operators.Chebyshev, v,
                            u.get_sampling_type(main_cfg["sampling"]), WR.check_cs)

	### Execute test
    test_result = test(spde_model, wr_model, main_cfg, sparse_cfg, CrossCheck(main_cfg["nb_tests"]), test_result)
    # test_result = test(spde_model, wr_model, main_cfg, [CrossCheck(main_cfg["nb_tests"])], main_cfg["experiment_name"], *test_result)


### Main
if __name__ == "__main__":
    
    parser = MLCSPGParse.define_MLCSPG_cli_parser()
    args = parser.parse_args(sys.argv[1:])
    main_cfg, pde_cfg, sparse_cfg = MLCSPGParse.parse_cfg_file_and_args(args)
	
    
    Main(main_cfg, pde_cfg, sparse_cfg)