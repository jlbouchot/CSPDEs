import argparse
import configparser
from pathlib import Path
from typing import Any, Dict, Tuple

int_arg_list = ["nb_level", "nb_iter", "nb_tests", "sJ", "sL", "nb_cosines", "mesh_x", "mesh_y", "l_start", "ansatz_space"]
float_arg_list = ["const_sJ", "exponent", "abar", "fluctuation_importance", "gamma", "tol_res", "power", "weight_cosine", 
        "dat_constant", "t_0", "t_prime", "smooth_0", "smooth_t"]
bool_arg_list = ["better_compute"]
args_list = ["output_file", "nb_cosines", "mesh_x", "mesh_y", "nb_level", "recovery_algo", "gamma", "l_start", 
        "sampling", "nb_iter", "tol_res", "nb_tests", "power", "abar", "fluctuation_importance", "weight_cosine", 
        "dat_constant", "prefix_precompute", "better_compute", "ansatz_space", "t_0", "t_prime", "smooth_0", "smooth_t", 
        "const_sJ", "no_compute", "exponent"]


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
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputDiffusionMLPolynomial", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-x", "--mesh-x", help="Size of the coarsest level (number of grid points) in the x direction", default=2000, required=False)
    parser.add_argument("-y", "--mesh-y", help="Size of the coarsest level (number of grid points) in the y direction", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=50, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default="pragmatic", required=False)
    parser.add_argument("-n", "--nb-tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-p", "--power", help="Power of the decay of the trigonometric expansion (~ mu)", default=4.0, required=False)
    parser.add_argument("-a", "--abar", help="Value of the mean field", default=10, required=False)
    parser.add_argument("-b", "--better-compute", help="Should the computations be done on the fly, using tensor representation (Default is TRUE)", default="True", required=False)
    parser.add_argument("-c", "--dat_constant", help="Multiplicative constant for expression of s_L", default=15., required=False)
    parser.add_argument("-d", "--nb-cosines", help="Number of random cosine and sine parameters", default=5, required=False)
    parser.add_argument("-e", "--tol-res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=1e-4, required=False)
    parser.add_argument("-E", "--exponent", help="Power of the polynomial weight (~ alpha)", default=1.0/4.0, required=False)
    parser.add_argument("-f", "--prefix-precompute", help="How should the precomputed data for this test be called?", default="testingWCosine", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant coefficients", default=1.035, required=False)
    parser.add_argument("-i", "--fluctuation-importance", help="What is the importance of the fluctuations with respect to the mean field (default is 1)", default=1, required=False)
    parser.add_argument("-j", "--ansatz-space", help="What type of Ansatz space is used? (Default is 0)", default="0", required=False)
    parser.add_argument("-k", "--no-compute", help="Should we skip all computations and only check values for s, m, and N (Default = False)", default=False, required=False)
    parser.add_argument("-w", "--weight-cosine", help="How much weight the local cosine carries (Default = 1, ~ Upsilon)", default=1, required=False)
    parser.add_argument("--t_0", help="What is the smoothness of the data (Default is 1)", default="1", required=False)
    parser.add_argument("--t_prime", help="What is the smoothness of the functional (Default is 1)", default="1", required=False)
    parser.add_argument("--smooth_0", help="What kind of smoothness in the original space can be expected (Default is 1/2)", default="0.5", required=False)
    parser.add_argument("--smooth_t", help="What kind of smoothness in the smooth space can be expected (Default is 1/2)", default="0.5", required=False)
    parser.add_argument("--const_sJ", help="What is the expected constant in the expression of s_J (Default is 15)", default="15", required=False)
    # Add config file parameter
    return parser

def parse_cli_args(config, cli_args, field_list): 
    '''
        TODO: Add the other config things (pde and sparse solvers)
    '''
    for field in field_list:
        '''
        [
            "input_mesh", "output_file",
            "max_iter", "tolerance",
            "levels", "alpha", "beta"
        ]:
        '''
        val = getattr(cli_args, field, None)
        if val is not None:
            config[field] = val

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
    parse_config_file(main_cfg, pde_solver_cfg, sparse_solver_cfg, path_to_file = "mlcspg-default-cfg.ini")
    # Check if a specific config file exists for this particular experiment
    specific_fname = getattr(all_args, "cfg", None)
    if specific_fname:
        parse_config_file(main_cfg, pde_solver_cfg, sparse_solver_cfg, path_to_file = specific_fname)
    # Parse all CLI argument and overload previous values if conflicting
    parse_cli_args(main_cfg, all_args, args_list) # TODO: Add the other config dicts
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
