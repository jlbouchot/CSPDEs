import WR

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML2D_tiny, PiecewiseConstantDiffusionFEMModelML2D_small, PiecewiseConstantDiffusionFEMModelML2D, PiecewiseConstantDiffusionFEMModelML2D_large, TrigCoefficient, ConstantCoefficient, Average

from Check_ML import test, CrossCheck

import sys
import numpy as np


def Main(outfile):
    
    ## SPDEModel
    d     = 9

    ## Parameters -- Have to be parts of the command line to be done neatly
    variability = 2.0
    grid_points = tuple([100,100])
    abar = 5.0
    nb_iter = 50
    epsilon     = 1e-8
    L = 3	
    num_tests = 20
    dat_constant = 8
    gamma = 1.1 # Used for the uniform weights

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    d = 9
    spde_model = PiecewiseConstantDiffusionFEMModelML2D_tiny(abar, ConstantCoefficient(1.0), variability, grid_points, Average())

    test_result_tiny = outfile + '_2D_tiny', None
    v = np.hstack((np.repeat(gamma, d), [np.inf]))
    wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)
    test(spde_model, wr_model, nb_iter, epsilon, L, [CrossCheck(num_tests)], dat_constant, "piecewiseLinearConstant2D", test_result_tiny)


    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    d = 16
    spde_model = PiecewiseConstantDiffusionFEMModelML2D_small(abar, ConstantCoefficient(1.0), variability, grid_points, Average())

    test_result_small = outfile + '_2D_small', None
    v = np.hstack((np.repeat(gamma, d), [np.inf]))
    wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)
    test(spde_model, wr_model, nb_iter, epsilon, L, [CrossCheck(num_tests)], dat_constant, test_result_small)

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    d = 25
    spde_model = PiecewiseConstantDiffusionFEMModelML2D(abar, ConstantCoefficient(1.0), variability, grid_points, Average())

    test_result = outfile + '_2D', None
    v = np.hstack((np.repeat(gamma, d), [np.inf]))
    wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)
    test(spde_model, wr_model, nb_iter, epsilon, L, [CrossCheck(num_tests)], dat_constant, test_result)

    # Create FEMModel with given diffusion coefficient, goal functional and initial mesh size
    d = 36
    spde_model = PiecewiseConstantDiffusionFEMModelML2D_large(abar, ConstantCoefficient(1.0), variability, grid_points, Average())

    test_result_large = outfile + '_2D_large', None
    v = np.hstack((np.repeat(gamma, d), [np.inf]))
    wr_model   = WR.WRModel(WR.Algorithms.whtp, WR.Operators.Chebyshev, v,
                                    WR.cs_pragmatic_m, WR.check_cs)
    test(spde_model, wr_model, nb_iter, epsilon, L, [CrossCheck(num_tests)], dat_constant, test_result_large)

            

### Main
if __name__ == "__main__":
    Main(sys.argv[1])
