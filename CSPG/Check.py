import numpy as np
import datetime
import time
import shelve

from collections import namedtuple

from CSPDE import CSPDE

TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 's', 'cspde_result'])
def test(spde_model, wr_model, epsilon, s, checks = None, filename = None, cspde_result = None):
    ### Execute CSPDE algorithm
    cspde_result = CSPDE(spde_model, wr_model, epsilon, s, cspde_result)


    ### Output results and check
    print("\nPostprocessing and outputting solution ...")

    ## Save results
    dt = datetime.datetime.fromtimestamp(time.time()).isoformat()
    if filename is None:
        filename = 'results_{0}'.format(dt)

    print("   Writing results to {0} ...".format(filename))
    d     = shelve.open(filename)
    d[dt] = TestResult(spde_model, wr_model, epsilon, s, cspde_result)
    d.close()

    ## Execute checks
    if not checks is None:
        print("   Executing checks ... ")
        map(lambda C: C(spde_model, wr_model, epsilon, cspde_result), checks)


    return filename, cspde_result


class CrossCheck:
    def __init__(self, num_tests):
        self.num_tests = num_tests

    def __call__(self, spde_model, wr_model, epsilon, cspde_result, y_truth=None):
        # Compute truth values
        Z_cross = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (self.num_tests, cspde_result.d)))

        if y_truth is None:
            # Compute truth values of functional at new samples
            y_truth = spde_model.samples(Z_cross, epsilon/2.)

        # Compute reconstructed values
        F_recon = wr_model.operator.create(cspde_result.J_s, Z_cross, normalization=np.sqrt(cspde_result.m))
        y_recon = F_recon.apply(cspde_result.result.x)

        # Compare
        difference = np.abs(y_truth - y_recon)
        print("Maximum error: {0}; Average error: {1}".format(difference.max(), difference.sum()/self.num_tests))

        return y_truth
