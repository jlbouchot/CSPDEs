import numpy as np
import datetime
import time
import shelve

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@bit.edu.cn"
__status__ = "Development"
__lastmodified__ = "2019/03/06"

from collections import namedtuple

from CSPDE_ML import CSPDE_ML

TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
def test(spde_model, wr_model, nb_iter, epsilon, J, L, checks = None, dat_constant = 10, ansatz_space = 0, prefix_fname = None, filename = None, cspde_result = None):
    ### Execute CSPDE algorithm
    sampling_fname = prefix_fname + 'sampling_points_Lmax' + str(L)
    datamtx_fname = prefix_fname + 'datamtx_Lmax' + str(L)
    

    cspde_result = CSPDE_ML(spde_model, wr_model, nb_iter, epsilon, J, L, dat_constant, ansatz_space, cspde_result, sampling_fname, datamtx_fname)

    ### Output results and check
    print("\nPostprocessing and outputting solution ...")

    ## Save results
    dt = datetime.datetime.fromtimestamp(time.time()).isoformat()
    if filename is None:
        filename = 'results_{0}'.format(dt)

    print("   Writing results to {0} ...".format(filename))
    d     = shelve.open(filename)
    d[dt] = TestResult(spde_model, wr_model, epsilon, L, cspde_result)
    d.close()

    ## Execute checks
    if not checks is None:
        print("   Executing checks ... ")
        list(map(lambda C: C(spde_model, wr_model, nb_iter, epsilon, cspde_result), checks)) # the "list(...)" is required due to the new Python 3.x updates

    ## Compute (L_2 estimation of) true Cheb coefs. 

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
