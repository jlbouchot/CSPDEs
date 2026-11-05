'''
Test file for running a PDE on a "2 dimensional grid" using the new dictionary propagation of information. 
The old tests will be legacy, deprecated, and then completely removed eventually.

This new test file tests the split of the various modeling elements of the PDE model:
* a discretization class, handling all sorts of mesh related issues (and function spaces and boundary conditions). TODO: bc could be potentially handled elsewhere if we decide to extend the approach to more complex bc.
* a physics class, handling the PDE and its representation in the weak form, and the goal functional.
* a solver class, adapted to various problems (linear / nonlinear solvers / assembled or not / ...)
* a sampler class, responsible for generating the samples of the solution, given a set of parameters. It makes the link between all the previous classes. 

The problem handles an affine linear parametric diffusion problem with CS recovery using polynomial weights. 
'''

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/29"
__lastmodified__ = "2026/11/04"

import os.path
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

PDE_PRESET_VALUES = (("OMP_NUM_THREADS", "omp_threads"),
                     ("OPENBLAS_NUM_THREADS", "blas_threads"),
                     ("MKL_NUM_THREADS", "blas_threads"))
SPARSE_PRESET_VALUES = (("NUMBA_NUM_THREADS", "numba_threads"),)

# --- Configuration first: thread counts must be in the environment before numpy / numba / dolfin are loaded ---
import utilities.MLCSPG_param_parser as MLCSPGParse

if __name__ == "__main__":
    _args = MLCSPGParse.define_MLCSPG_cli_parser().parse_args(sys.argv[1:])
    MAIN_CFG, PDE_CFG, SPARSE_CFG = MLCSPGParse.parse_cfg_file_and_args(_args)
    for var, key in PDE_PRESET_VALUES:
        os.environ[var] = str(PDE_CFG.get(key, 1))     # str: environ only takes strings; plain assignment: the ini wins
    for var, key in SPARSE_PRESET_VALUES:
        os.environ[var] = str(SPARSE_CFG.get(key, 1))

# --- Only now the heavy imports ---
import WR
import utilities.utils as u
from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModelML, WeightedCosine2D, ConstantCoefficient, Average, AffineExpansionFEMModelML, WeightedCosSin2D
from Check_ML import test, CrossCheck
import numpy as np
from SPDESampler import Physics
from SPDESampler import Samplers
from SPDESampler import Discretisation
from SPDESampler import Functionals
from SPDESampler.Physics import Diffusion
from SPDESampler.Samplers import AffineLinearSampler
from SPDESampler.Discretisation import UnitDiscretisation
from SPDESampler.Functionals import AverageQOI

def Main(main_cfg, pde_cfg, sparse_cfg): 

    # Define the spde model and all its physics / solver / etc ...
    
    # grid_points = tuple(int(2**(main_cfg["l_start"])*dummy) for dummy in [main_cfg["mesh_x"], main_cfg["mesh_y"], main_cfg["mesh_z"]][0:main_cfg["n"]])
    mesh_size = tuple(int(dummy) for dummy in [main_cfg["mesh_x"], main_cfg["mesh_y"], main_cfg["mesh_z"]][0:main_cfg["n"]])
    		
    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    do_cossin = True
    if do_cossin:
        coef_expansion = WeightedCosSin2D(main_cfg["nb_cosines"]*2, main_cfg["power"], main_cfg["fluctuation_importance"], main_cfg["abar"], main_cfg["weight_cosine"])
    else:
        coef_expansion = WeightedCosine2D(main_cfg["nb_cosines"], main_cfg["power"], main_cfg["fluctuation_importance"], main_cfg["abar"], main_cfg["weight_cosine"])
    #spde_model = DiffusionFEMModelML(coef_expansion, ConstantCoefficient(10.0),
    #                                   Average(), grid_points, pde_cfg) 
    # spde_model = AffineExpansionFEMModelML(coef_expansion, ConstantCoefficient(10.0),
    #                                    Average(), grid_points, pde_cfg) 
    phys = Diffusion(coef_expansion)
    # TODO: 'cholmod' should be replaced with some informations gathered from the ini files via the pde_cfg.
    # TODO: somehow, we should be careful with handling parameters of the numerical solver. Some iterative solvers may have parameters needed later on. Should we initialise a solver here? Or should we pass the parameters downwards and initialise at first call?
    spde_model = [phys, AffineLinearSampler, UnitDiscretisation, AverageQOI(), pde_cfg["linear_solver"], mesh_size]
                                       
    ### Reconstruction Model
    # TODO: Rename things correctly and pass them as parameters
    v = np.hstack((main_cfg["gamma"]*np.power([val+1 for val in range(main_cfg["nb_cosines"]) for _ in ((0,1) if do_cossin else (0,)) ], main_cfg["exponent"]), [np.inf]))

    operator = u.get_operator_type(main_cfg["operator"])
    
    if not u.is_tensor_based_possible(sparse_cfg["recovery_algo"]):
       main_cfg["operator"] = "Chebyshev"
       operator =  WR.Operators.Chebyshev
    # TODO: Remove the 'if' from this test file towards the more "mechanical" parts of the code.
    wr_model = WR.WRModel(sparse_cfg["recovery_algo"], operator, v, 
                        u.get_sampling_type(main_cfg["sampling"]), WR.check_cs)

	### Execute test
    test_result = test(spde_model, wr_model, main_cfg, sparse_cfg, pde_cfg, CrossCheck(main_cfg["nb_tests"]), True)

### Main
if __name__ == "__main__":
    Main(MAIN_CFG, PDE_CFG, SPARSE_CFG)
