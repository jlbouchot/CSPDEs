import MLCSPG_param_parser as MLCSPGParse

import sys
import argparse


__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019-2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/01/16"



# def Main(outfile, d = 10, L_max = 4, orig_mesh_size = 2000):
def Main(main_cfg, pde_cfg, sparse_cfg):

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

### Main
if __name__ == "__main__":
    
    parser = MLCSPGParse.define_MLCSPG_cli_parser()
    args = parser.parse_args(sys.argv[1:])
    main_cfg, pde_cfg, sparse_cfg = MLCSPGParse.parse_cfg_file_and_args(args)
	
    
    Main(main_cfg, pde_cfg, sparse_cfg)
    # Main(sys.argv[1])
