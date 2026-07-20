'''
Test file for running a PDE on a "1 dimensional grid" using the new dictionary propagation of information. 
The old test_wCosineXXX will be legacy, deprecated, and then completely removed eventually.

The problem handles a parametric diffusion problem with CS recovery using polynomial weights. 
'''

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/07/09"
__lastmodified__ = "2026/07/17"

import os.path
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import utilities.MLCSPG_param_parser as MLCSPGParse

import WR

import utilities.utils as u

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, WeightedCosine1D, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import numpy as np

def Main(main_cfg, pde_cfg, sparse_cfg): 
    
    grid_points = tuple(int(2**(main_cfg["l_start"])*dummy) for dummy in [main_cfg["mesh_x"], main_cfg["mesh_y"], main_cfg["mesh_z"]][0:main_cfg["n"]])
    		
    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    coef_expansion = WeightedCosine1D(main_cfg["nb_cosines"], main_cfg["power"], main_cfg["fluctuation_importance"], main_cfg["abar"], main_cfg["weight_cosine"])
    spde_model = DiffusionFEMModelML(coef_expansion, ConstantCoefficient(10.0),
                                       Average(), grid_points, pde_cfg) 
                                       
    test_result = None
    ### Reconstruction Model
    # TODO: Rename things correctly and pass them as parameters
    v = np.hstack((main_cfg["gamma"]*np.power([val+1 for val in range(main_cfg["nb_cosines"]) for dummy_variable in (0,1)], main_cfg["exponent"]), [np.inf]))

    if not u.is_tensor_based_possible(sparse_cfg["recovery_algo"]):
       main_cfg["do_tensor"] = False 
    # TODO: Remove the 'if' from this test file towards the more "mechanical" parts of the code.
    if main_cfg["do_tensor"]: 
        wr_model   = WR.WRModel(sparse_cfg["recovery_algo"], WR.Operators.Cheb_Alt, v, 
                            u.get_sampling_type(main_cfg["sampling"]), WR.check_cs)
    else: # The basic way.
        wr_model   = WR.WRModel(sparse_cfg["recovery_algo"], WR.Operators.Chebyshev, v,
                            u.get_sampling_type(main_cfg["sampling"]), WR.check_cs)

	### Execute test
    test_result = test(spde_model, wr_model, main_cfg, sparse_cfg, pde_cfg, CrossCheck(main_cfg["nb_tests"]), test_result)

### Main
if __name__ == "__main__":
    parser = MLCSPGParse.define_MLCSPG_cli_parser()
    args = parser.parse_args(sys.argv[1:])
    main_cfg, pde_cfg, sparse_cfg = MLCSPGParse.parse_cfg_file_and_args(args)
	
    Main(main_cfg, pde_cfg, sparse_cfg)
