import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML, LinearCoefficient, ConstantCoefficient, Integration

from Check import test, CrossCheck

import sys
import numpy as np


def Main(outfile):

    ## SPDEModel
    d  = 9

    mesh_size = d*200 # corresponds to 1800 grid points for our example

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    abar = 5.0
    variability = 1.0/(d+1)
    upper_bound = abar + variability
    a  = [LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
          LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound),
          LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound), LinearCoefficient(abar, upper_bound)] # Make sure you have d of those
    spde_model = PiecewiseConstantDiffusionFEMModelML(a, ConstantCoefficient(1.0), mesh_size)

    test_result = outfile, None

    for s in range(35,156,15):
        for gamma in np.arange(1.03, 1.08, 0.01)[::-1]:
            ## Reconstruction Model
            v = [gamma, gamma, gamma, gamma, gamma, np.inf]

            wr_model = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                  WR.cs_pragmatic_m, WR.check_cs)

            ### Execute test
            test_result = test(spde_model, wr_model, epsilon, s, [], *test_result)


### Main
if __name__ == "__main__":
    Main(sys.argv[1])
