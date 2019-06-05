import numpy as np
import datetime
import time
import shelve
import os
import json # helps keep the config file easily readable (and tunable)!

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@bit.edu.cn"
__status__ = "Development"
__lastmodified__ = "2019/06/05"

from collections import namedtuple

from CSPDE_ML import CSPDE_ML

TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
# def test(spde_model, wr_model, nb_iter, epsilon, J, L, checks = None, dat_constant = 10, p = 2/3, p0 = 1/3, t = 1, tprime = 1, energy_constant = 10, ansatz_space = 0, prefix_fname = None, filename = None, cspde_result = None):
def test(spde_model, wr_model, dict_config, checks = None, prefix_fname = None, filename = None, cspde_result = None):

    nb_iter = dict_config["iter"] 
    epsilon = dict_config["tolres"]
    J = dict_config["J"] 
    L = dict_config["L"] 
    dat_constant = dict_config["s_L"] 
    p = dict_config["p"]
    p0 = dict_config["p0"] 
    t = dict_config["t"]
    tprime = dict_config["tprime"]
    energy_constant = dict_config["s_J"]
    ansatz_space = dict_config["ansatz"] 


    # Create target Directory if don't exist
    if not os.path.exists(prefix_fname):
        os.mkdir(prefix_fname)
    else:    
        print("WARNING: Directory " , prefix_fname ,  " already exists \n ---> You might overwrite important files!!")

    with open(os.path.join(prefix_fname, "config_file.txt"),'w') as f_handler:
        json.dump(dict_config,f_handler)

    ### Execute CSPDE algorithm
    sampling_fname = os.path.join(prefix_fname, 'sampling_points_Lmax' + str(L))
    datamtx_fname = os.path.join(prefix_fname, 'datamtx_Lmax' + str(L))
    

    cspde_result = CSPDE_ML(spde_model, wr_model, dict_config, cspde_result, sampling_fname, datamtx_fname)
#     cspde_result = CSPDE_ML(spde_model, wr_model, nb_iter, epsilon, J, L, dat_constant, ansatz_space, cspde_result, sampling_fname, datamtx_fname, t, tprime, p0, p, energy_constant)

    ### Output results and check
    print("\nPostprocessing and outputting solution ...")

    ## Save results
    dt = datetime.datetime.fromtimestamp(time.time()).isoformat()
    if filename is None:
        filename = 'results_{0}'.format(dt)

    print("   Writing results to {0} ...".format(filename))
    d     = shelve.open(os.path.join(prefix_fname,filename))
    d[dt] = TestResult(spde_model, wr_model, epsilon, L, cspde_result)
    d.close()

    ## Execute checks
    if not checks is None:
        print("   Executing checks ... ")
        list(map(lambda C: C(spde_model, wr_model, nb_iter, epsilon, cspde_result), checks)) # the "list(...)" is required due to the new Python 3.x updates


    return filename, cspde_result

class CrossCheck:
    def __init__(self, num_tests):
        self.num_tests = num_tests

    def __call__(self, spde_model, wr_model, nb_iter, epsilon, cspde_result, y_truth=None):
        # Compute truth values
        if self.num_tests == [] or self.num_tests is None: 
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
