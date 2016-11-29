from plot_comparison_1D_ML import plot_1D
from plot_comparison_1D_ML import plot_polynomialweights_1D as ppw1D
from plot_convergence_hL import plot_1D_hL
from plot_coefValues import plot_coef_values as pcv

import numpy as np


## First plot the convergence for the small and large variances of the piecewise constant diffusion case in 1D
# Make sure we have the appropriate file names
file_small_var = 'pwldML_whtp_var_small_L3_grid200_g105'
file_big_var = 'pwldML_whtp_var_big_L3_grid200_g105'
# Define the important variables we want to display
v_to_show = [1.05]
param_to_calculate = [0,2,4,6]
L_to_show = [1,2,3]
# Let's go with the graphs
plot_1D(file_small_var, 'plot_pointwiseConvergence_pwld1D_smallVar.png', L_to_show, v_to_show, param_to_calculate)
plot_1D(file_big_var, 'plot_pointwiseConvergence_pwld1D_bigVar.png', L_to_show, v_to_show, param_to_calculate)

### Same thing with the 2D piecewise constant diffusion coefficients  
## Make sure we have the appropriate file names
#file_small_var = 'pwld2D_grid200_var_small_L3_g107' 
#file_big_var = 'pwld2D_grid200_var_big_L3_g107' 
## Define the important variables we want to display
#v_to_show = [1.07] 
#param_to_calculate = [0,6,12,18] 
#L_to_show = [1,2,3] 
## Let's go with the graphs
#plot_1D(file_small_var, 'plot_pointwiseConvergence_pwld2D_smallVar.png', L_to_show, v_to_show, param_to_calculate) 
#plot_1D(file_big_var, 'plot_pointwiseConvergence_pwld2D_bigVar.png', L_to_show, v_to_show, param_to_calculate) 


## Pointwise convergence for a trigonometric expansion of the diffusion coefficient
# Define the important files
from plot_comparison_1D_ML import plot_polynomialweights_1D as ppw1D
from plot_convergence_hL import plot_1D_hL

file_polynomial = 'diffML_trig_bp_poly_d40_grid2000_L4_datconstant8_c102' 
# What to display
param_to_calculate = [0,2,4,6] 
L_to_show = [2,3,4] 
# Let's go with the graphs
# ppw1D(file_polynomial, 'plot_pointwiseConvergence_Polynomial.png', L_to_show, param_to_calculate) 
plot_1D_hL('diffML_trig_bp_poly_d40_datconstant8_c102_fixedhL20000', 'plot_fixhL20000_c102_const8_bp_l234.png', L_to_show, param_to_calculate)


## Plot target accuracy vs xnumber of levels for the trigonometric expansion
# Define the generic filename
gen_fname = 'diffML_whtp_d5_g1055p_40000' # Check this! (and then delete comment)
v_to_show = [1.055] # Check this! (and then delete comment)
param_to_calculate = [0,1,2,3] # Check this! (and then delete comment)
L_to_show = [1,2,3] # Check this! (and then delete comment)
plot_1D_hL('diffML_whtp_d5_g1055p_40000', 'plot_fixhL.png', L_to_show, v_to_show, param_to_calculate) # Change this accordingly (and delete once done)

## Compare true coefs with estimated coefs
# piecewise constant diffusion coefficients in 2D
pwld2D_file = 'pwld2D_grid200_var_big_L3' 
figname = 'comp_true_coefs_pwld2D_big_L3.png'


from plot_comparison_1D_ML import plot_1D
from plot_convergence_hL import plot_1D_hL
from plot_coefValues import plot_coef_values as pcv

import numpy as np

pwld2D_file = 'smallTestsPoisson' 
figname = 'coefPoissonTests.png'
coefs_to_display_ratio = 0.1 # hopefully, this should be more than enough
pcv(pwld2D_file, figname, coefs_to_display_ratio)

diffML_file = 'diffML_trig_whtp_d30_g11_grid2000_L3_large'
figname = 'comp_true_coefs_diffML_trig_g11_d30_grid200_big_L3.png'
pcv(diffML_file, figname, coefs_to_display_ratio)


## And another example
from plot_comparison_1D_ML import plot_1D
from plot_convergence_hL import plot_1D_hL
from plot_coefValues import plot_coef_values as pcv

import numpy as np

pwld2D_file = 'smallDim6_PWLD_whtp_var_small_grid200_L4_g105' 
figname = 'smallDim6_pwld_compareCoefs.png'
coefs_to_display_ratio = 0.08 # hopefully, this should be more than enough
pcv(pwld2D_file, figname, coefs_to_display_ratio)


# plot_1D('pwldML_whtp_var_1', 'plot_pwld.png', L_to_show, v_to_show, param_to_calculate)


from collections import namedtuple

MCResults = namedtuple('MCResults', ['nu', 'MC_est', 'l1_est'])
from plot_coefValues import plot_coef_from_file as pcff

datafile = 'smallDim6_PWLD_whtp_var_small_grid2000_L3_g105_datconstant20_MCResults'
figname = 'compare_true_coefs.png'
pcff(datafile,figname)
