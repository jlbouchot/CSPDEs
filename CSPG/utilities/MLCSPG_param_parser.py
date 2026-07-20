import argparse
import configparser
from pathlib import Path
import os.path
import inspect
from typing import Any, Dict, Tuple

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/07/09"
__lastmodified__ = "2026/07/09"

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)


int_arg_list = ["nb_level", "nb_iter", "nb_tests", "sJ", "sL", "nb_cosines", "mesh_x", "mesh_y", "l_start", "ansatz_space", "n", "degree"]
float_arg_list = ["const_sj", "exponent", "abar", "fluctuation_importance", "gamma", "tol_res", "power", "weight_cosine", 
        "dat_constant", "t_0", "t_prime", "p_0", "p_t"]
bool_arg_list = ["do_tensor", "no_compute"]
args_list = ["output_file", "nb_cosines", "mesh_x", "mesh_y", "nb_level", "recovery_algo", "gamma", "l_start", 
        "sampling", "nb_iter", "tol_res", "nb_tests", "power", "abar", "fluctuation_importance", "weight_cosine", 
        "dat_constant", "prefix_precompute", "do_tensor", "ansatz_space", "t_0", "t_prime", "p_0", "p_t", 
        "const_sj", "no_compute", "exponent", "preconditioner", "linear_solver", "n", "degree", "experiment_name", "elements"]
pde_arg_list = ["preconditioner", "linear_solver", "elements", "degree"]
sparse_arg_list = ["recovery_algo", "nb_iter", "tol_res"]
main_arg_list = ["output_file", "nb_cosines", "mesh_x", "mesh_y", "nb_level", "l_start", "sampling", "power", "abar", 
        "fluctuation_importance", "weight_cosine", "dat_constant", "prefix_precompute", "do_tensor", "ansatz_space", 
        "t_0", "t_prime", "p_0", "p_t", "const_sj", "no_compute", "exponent", "n", "experiment_name"]


# def parse_config_file(path_to_file = None: str, exp = None: str) -> Tuple[Dict[str,Any], Dict[str,Any]]:
def parse_config_file(main_cfg, pde_solver_cfg, sparse_solver_cfg, path_to_file = None):
    if path_to_file:
        parser = configparser.ConfigParser()
        config_path = Path(path_to_file)
        if config_path.is_file():
            parser.read(config_path)
        else:
            raise FileNotFoundError(f"Config file {config_path} not found")

    if parser:
        if "main" in parser:
            for key, val in parser["main"].items():
                main_cfg[key] = val
        
        if "pdesolver" in parser:
            for key, val in parser["pdesolver"].items():
                pde_solver_cfg[key] = val

        if "sparsesolver" in parser: 
            for key, val in parser["sparsesolver"].items():
                sparse_solver_cfg[key] = val


def define_MLCSPG_cli_parser(): 
    ''' List all potential command line arguments to be used in the MLCSPG algorithm. '''
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-o", "--output_file", help="File to write the results", default=None, required=False)
    parser.add_argument("-L", "--nb_level", help="Number of levels used", default=None, required=False)
    parser.add_argument("-x", "--mesh_x", help="Size of the coarsest level (number of grid points) in the x direction", default=None, required=False)
    parser.add_argument("-y", "--mesh_y", help="Size of the coarsest level (number of grid points) in the y direction", default=None, required=False)
    parser.add_argument("-N", "--nb_iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=None, required=False)
    parser.add_argument("-r", "--recovery_algo", help="String for the algorithm for weighted l1 recovery", default=None, required=False)
    parser.add_argument("-s", "--l_start", help="Instead of going through all the levels, give it a starting point", default=None, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default=None, required=False)
    parser.add_argument("-n", "--nb_tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-p", "--power", help="Power of the decay of the trigonometric expansion (~ mu)", default=None, required=False)
    parser.add_argument("-a", "--abar", help="Value of the mean field", default=None, required=False)
    parser.add_argument("-b", "--do_tensor", help="Should the computations be done on the fly, using tensor representation (Default is TRUE)", default=None, required=False)
    parser.add_argument("-c", "--dat_constant", help="Multiplicative constant for expression of s_L", default=None, required=False)
    parser.add_argument("-d", "--nb_cosines", help="Number of random cosine and sine parameters", default=None, required=False)
    parser.add_argument("-e", "--tol_res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=None, required=False)
    parser.add_argument("-E", "--exponent", help="Power of the polynomial weight (~ alpha)", default=None, required=False)
    parser.add_argument("-f", "--experiment_name", help="How should the precomputed data for this test be called?", default=None, required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant coefficients", default=None, required=False)
    parser.add_argument("-i", "--fluctuation_importance", help="What is the importance of the fluctuations with respect to the mean field (default is 1)", default=None, required=False)
    parser.add_argument("-j", "--ansatz_space", help="What type of Ansatz space is used? (Default is 0)", default=None, required=False)
    parser.add_argument("-k", "--no_compute", help="Should we skip all computations and only check values for s, m, and N (Default = False)", default=None, required=False)
    parser.add_argument("-w", "--weight_cosine", help="How much weight the local cosine carries (Default = 1, ~ Upsilon)", default=None, required=False)
    parser.add_argument("--t_0", help="What is the smoothness of the data (Default is 1)", default=None, required=False)
    parser.add_argument("--t_prime", help="What is the smoothness of the functional (Default is 1)", default=None, required=False)
    parser.add_argument("--p_0", help="What kind of smoothness in the original space can be expected (Default is 1/2)", default=None, required=False)
    parser.add_argument("--p_t", help="What kind of smoothness in the smooth space can be expected (Default is 1/2)", default=None, required=False)
    parser.add_argument("--const_sJ", help="What is the expected constant in the expression of s_J (Default is 15)", default=None, required=False)
    parser.add_argument("--preconditioner", help="Preconditioner chosen among those available from FEniCS", default=None, required=False) 
    parser.add_argument("--linear_solver", help="Type of solver used for the PDE solves", default=None, required=False)
    parser.add_argument("--cfg", help="Specific config file for the current experiment", default=None,required=False)
    parser.add_argument("--degree", help="Degree of the finite elements used in the FEM discretization", default=1, required=False)
    parser.add_argument("--elements", help="Type of finite elements used in the FEM discretization", default="Lagrange", required=False)
    # Add config file parameter
    return parser

def parse_cli_args(main_cfg, pde_solver_cfg, sparse_solver_cfg, cli_args): 
    '''
        TODO: Potentially add the possibility to have more than one config file, i.e. iterate over specific names
    '''
    # Check if a specific config file exists for this particular experiment
    specific_fname = getattr(cli_args, "cfg", None)
    if specific_fname:
        parse_config_file(main_cfg, pde_solver_cfg, sparse_solver_cfg, path_to_file = specific_fname)
    for field in main_arg_list:
        val = getattr(cli_args, field, None)
        if val is not None:
            main_cfg[field] = val
    for field in pde_arg_list:
        val = getattr(cli_args, field, None)
        if val is not None:
            pde_solver_cfg[field] = val
    for field in sparse_arg_list:
        val = getattr(cli_args, field, None)
        if val is not None:
            sparse_solver_cfg[field] = val

def parse_cfg_file_and_args(all_args): 
    '''
    all_args contains the namespace obtained from parsing the argv via the cli parser
    1. Parse the default config file
    2. Check if we have a specific config file within the command lines and overload the dicts if we have one
    3. Overload everything with the CLI arguments
    4. Cast numerics whenever needed
    '''
    main_cfg = {}
    pde_solver_cfg = {}
    sparse_solver_cfg = {}
    # Parse config file modifies in place main and solver cfg
    parse_config_file(main_cfg, pde_solver_cfg, sparse_solver_cfg, path_to_file = os.path.join(parentdir, 'data', 'mlcspg-default-cfg.ini'))
    # Parse all CLI argument and overload previous values if conflicting
    parse_cli_args(main_cfg, pde_solver_cfg, sparse_solver_cfg, all_args)
    main_keys = main_cfg.keys()
    pde_solver_keys = pde_solver_cfg.keys()
    sparse_solver_keys = sparse_solver_cfg.keys()
    for k in int_arg_list:
        if k in main_keys:
            main_cfg[k] = int(main_cfg[k])
        if k in pde_solver_keys:
            pde_solver_cfg[k] = int(pde_solver_cfg[k])
        if k in sparse_solver_keys:
            sparse_solver_cfg[k] = int(sparse_solver_cfg[k])

    for k in float_arg_list:
        if k in main_keys:
            main_cfg[k] = float(main_cfg[k])
        if k in pde_solver_keys:
            pde_solver_cfg[k] = float(pde_solver_cfg[k])
        if k in sparse_solver_keys:
            sparse_solver_cfg[k] = float(sparse_solver_cfg[k])
    
    for k in bool_arg_list:
        if k in main_keys:
            main_cfg[k] = main_cfg[k].lower() == "true"
        if k in pde_solver_keys:
            pde_solver_cfg[k] = pde_solver_cfg[k].lower() == "true"
        if k in sparse_solver_keys:
            sparse_solver_cfg[k] = sparse_solver_cfg[k].lower() == "true"

    return main_cfg, pde_solver_cfg, sparse_solver_cfg
