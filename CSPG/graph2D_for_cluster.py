# A little script for compraing our estimation with the ground truth on the cluster. 
# I have to do this way, because, for some reasons, it doesn't want to work... 
import numpy as np

from SPDE              import FEniCSModels
from SPDE.FEniCSModels import PiecewiseConstantDiffusionFEMModelML2D_small, LinearCoefficient, ConstantCoefficient, Integration, Average

step = 0.02 # Change back to 0.02
T = np.arange(-1+step, 1, step)
d  = 16
abar = 5.0
gamma = 1.1


variability = 1.0/25.0
spde_model = PiecewiseConstantDiffusionFEMModelML2D_small(abar, ConstantCoefficient(1.0), variability, tuple([1000,1000]), Average())

dump_fname = 'file_for_graph_2D_smallVar.npy'
dict_results_smallVar = {}
for k in [0,5,10,15]:
    Z = np.array([np.hstack((np.zeros(k), t, np.zeros(d-k-1))) for t in T])
    y_new = spde_model.samples(Z)
    dict_results_smallVar[k] = y_new
    np.save(dump_fname, dict_results_smallVar)


variability = 2.0
spde_model = PiecewiseConstantDiffusionFEMModelML2D_small(abar, ConstantCoefficient(1.0), variability, tuple([1000,1000]), Average())

dump_fname = 'file_for_graph_2D_largeVar.npy'
dict_results_bigVar = {}
for k in [0,5,10,15]:
    Z = np.array([np.hstack((np.zeros(k), t, np.zeros(d-k-1))) for t in T])
    y_new = spde_model.samples(Z)
    dict_results_bigVar[k] = y_new
    np.save(dump_fname, dict_results_smallVar)


