from plot_comparison_1D_ML import plot_1D
from plot_convergence_hL import plot_1D_hL

import numpy as np

# v_to_show = [1.055]
# param_to_calculate = [0,2,4,6]
# L_to_show = [1,2,3]

# plot_1D('diffML_whtp_d5_g1055p_2000', 'plot_fixh0.png', L_to_show, v_to_show, param_to_calculate)

v_to_show = [1.055]
param_to_calculate = [0,2,4,6]
L_to_show = [1,2,3,4]
plot_1D_hL('diffML_whtp_d5_g1055p_40000', 'plot_fixhL.png', L_to_show, v_to_show, param_to_calculate)

v_to_show = [1.05]
param_to_calculate = [0,1,2,3]
L_to_show = [1,2,3,4,5]
plot_1D('pwldML_whtp_var_1', 'plot_pwld.png', L_to_show, v_to_show, param_to_calculate)