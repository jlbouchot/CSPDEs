from plot_comparison_1D_ML import plot_1D

import numpy as np

v_to_show = [1.05]
param_to_calculate = [0,2,4,6]
L_to_show = [1,2,3,4]

plot_1D('pwldML_whtp_var_5', 'testGraphsPWLD.pdf', L_to_show, v_to_show, param_to_calculate)