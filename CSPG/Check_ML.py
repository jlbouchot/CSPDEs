import numpy as np
import datetime
import time
import shelve
import os
import json # helps keep the config file easily readable (and tunable)!

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015 - 2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/07/20"

from collections import namedtuple

from CSPDE_ML import CSPDE_ML

TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'L', 'cspde_result'])
def test(spde_model, wr_model, dict_config, sparse_config, pde_config, checks = None, cspde_result = None):

    J = dict_config["l_start"] 
    L = dict_config["nb_level"] 
    dat_constant = dict_config["dat_constant"]  # TODO: Eventually change this name! This is ugly. --> const_sL
    p = dict_config["p_t"]
    p0 = dict_config["p_0"] 
    t = dict_config["t_0"]
    tprime = dict_config["t_prime"]
    energy_constant = dict_config["const_sj"]
    ansatz_space = dict_config["ansatz_space"] 
    no_compute = dict_config["no_compute"]
    prefix_fname = dict_config["experiment_name"]
    filename = dict_config["output_file"]


    # Create target Directory if doesn't exist
    if not os.path.exists(prefix_fname):
        os.mkdir(prefix_fname)
    else:    
        print("WARNING: Directory " , prefix_fname ,  " already exists \n ---> You might overwrite important files!!")

    ### Execute CSPDE algorithm
    cspde_result = CSPDE_ML(spde_model, wr_model, dict_config, sparse_config, cspde_result)

    ### Output results and check
    print("\nPostprocessing and outputting solution ...")

    ## Save results
    dt = datetime.datetime.fromtimestamp(time.process_time()).isoformat()
    if filename is None:
        filename = 'results_{0}'.format(dt)

    if not no_compute:
        print("   Writing results to {0} ...".format(filename))
        d     = shelve.open(os.path.join(prefix_fname,filename))
        d[dt] = TestResult(spde_model, wr_model, L, cspde_result)
        d.close()

    with open(os.path.join(prefix_fname, filename + "_config_file.txt"),'w') as f_handler:
        f_handler.write("### MLCSPG configuration file ###\n")
        f_handler.write("##  MAIN DICTIONARY   ##\n")
        json.dump(dict_config,f_handler) # Probably no longer need this, since everything works with config files to start with.
        f_handler.write("\n\n")
        f_handler.write("##  SPARSE DICTIONARY ##\n")
        json.dump(sparse_config,f_handler)
        f_handler.write("\n\n")
        f_handler.write("##  PDE DICTIONARY    ##\n")
        json.dump(pde_config,f_handler) # Assuming you have a PDE config dictionary
        f_handler.write("\n")
    ## Execute checks
    if checks and (not no_compute):
        print("   Executing checks ... ")
        checks(spde_model, wr_model, cspde_result)
        # list(map(lambda C: C(spde_model, wr_model, nb_iter, epsilon, cspde_result), checks)) # the "list(...)" is required due to the new Python 3.x updates


    return filename, cspde_result

class CrossCheck:
    def __init__(self, num_tests):
        self.num_tests = num_tests

    def __call__(self, spde_model, wr_model, cspde_result, y_truth=None):
        # Compute truth values
        if self.num_tests == [] or self.num_tests is None or self.num_tests == 0: 
            return None
        else: 
            Z_cross = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (self.num_tests, cspde_result[0].d)))

            if y_truth is None:
                # Compute truth values of functional at new samples
                print("Computing true solutions")
                y_truth = spde_model.samples(Z_cross, 3) # consider a three times finer grid 

            # Compute reconstructed values
            y_recon = wr_model.estimate_ML_samples(cspde_result, Z_cross)
            spde_model.refine_mesh(0.3334)

            # Compare
            difference = np.abs(y_truth - y_recon)
            print("Maximum error: {0}; Average error: {1}".format(difference.max(), difference.sum()/self.num_tests))

            return y_truth
