import numpy as np
from correct_true_solution_PWLD import compute_true_avg as cta
from iterative_solution import compute_true_avg as iter_cta

from call_fenics_from_tests import compute_numerical_solutions as cns

# from dolfin import *


# import WR

# from SPDE              import FEniCSModels
# from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML, LinearCoefficient, ConstantCoefficient, Integration, Average

nb_ys = 100

# generate a random vector: 
# ys = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0],[0.5,0,0,0,0,0,0,0,0,0,0,0,0],[-0.5,0,0,0,0,0,0,0,0,0,0,0,0],[0.25,0,0,0,0,0,0,0,0,0,0,0,0],[0,0.5,0,0,0,0,0,0,0,0,0,0,0],[0,0,0.5,0,0,0,0,0,0,0,0,0,0]]).transpose()
# ys = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0]]).transpose()

# spde_model = PiecewiseConstantDiffusionFEMModelML(5, ConstantCoefficient(1.0), 1.0/14.0, 2000, Average())
# u = spde_model.solve(ys[:,0])
# 
# plot(u)
# interactive()
# 
# print assemble(u,Average())
# print "Chosen ys = \n", Z
# ys = np.random.uniform(-1,1, [nb_parameters,nb_ys])
# print "Former ys = \n", ys
for interesting_file in ['smallDim6_PWLD_whtp_var_small_grid2000_L4_g105', 'pwldML_whtp_var_small_L3_grid200_g105']:
# interesting_file = 'pwldML_whtp_var_small_L3_grid200_g105'
    numerical_solution, Z, finest_result = cns(interesting_file, nb_ys, None, 16)
    nb_parameters = Z.shape[0]

    rhs = 1
    variability = 1.0/(nb_parameters)
    mean_field = 5
    # avg_for_y = cta(ys, rhs, variability, mean_field)
    # print avg_for_y

    analytic_solution = iter_cta(Z.transpose(), rhs, variability, mean_field)

    difference = np.abs(numerical_solution - analytic_solution)
    print("Maximum error: {0}; Average error: {1} taken from {2} sample points".format(difference.max(), difference.sum()/nb_ys, nb_ys))

