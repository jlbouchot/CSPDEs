import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import DiffusionFEMModel, TrigCoefficient, ConstantCoefficient, Average

from Check import test, CrossCheck

import sys
import numpy as np


def Main(outfile):
    ### Input parameters
    epsilon     = 1e-4 # Change from 8 to test on quinoa

    ## SPDEModel
    d     = 2

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    spde_model = DiffusionFEMModel(TrigCoefficient(d, 1.0, 4.3), ConstantCoefficient(10.0),
                                       Average(), 100)

    test_result = outfile, None
    for s in range(35,156,15):
        for gamma in np.arange(1.03, 1.06, 0.01)[::-1]:
            ### Reconstruction Model
            v = np.hstack((np.repeat(gamma, 2*d), [np.inf]))

            wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)

            ## Number of tests
            num_tests = 2 # change from 10 for Quinoa tests

            ### Execute test
            test_result = test(spde_model, wr_model, epsilon, s, [CrossCheck(num_tests)], *test_result)


### Main
if __name__ == "__main__":
    Main(sys.argv[1])
