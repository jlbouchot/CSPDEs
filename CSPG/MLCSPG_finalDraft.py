from test_diff_avg_v_ML import Main as DIFF_V
from test_diff_avg_p_ML import Main as DIFF_P
from test_piecewiselineardiffusion_v import Main as PWLD_V # As in PieceWiseLinearDiffusion
from test_smallpiecewiselineardiffusion_v import Main as SPWLD
from test_dim5piecewiselineardiffusion_v import Main as PWLD5
from test_piecewiselineardiffusion_v_2D import Main as PWLD_2D
from test_piecewiselineardiffusion_v_2D_small import Main as PWLD_2D_small
from test_diff_avg_v_ML_2D import Main as DIFF_V_2D

from multiprocessing import Pool


# Compute solutions
p = Pool(4)

# p.apply_async(DIFF_V, ("diffML_trig_whtp_d20_g1055_grid2000", 10, tuple([2000]), 4, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(DIFF_V, ("diffML_trig_whtp_d30_g1055_grid2000", 15, tuple([2000]), 4, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(DIFF_V, ("diffML_trig_whtp_d30_g1055_fixedhL40000_L1", 15, tuple([40000]), 1, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(DIFF_V, ("diffML_trig_whtp_d30_g1055_fixedhL40000_L2", 15, tuple([20000]), 2, "whtp", 1.055, 2, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(DIFF_V, ("diffML_trig_whtp_d30_g1055_fixedhL40000_L3", 15, tuple([10000]), 3, "whtp", 1.055, 3, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(DIFF_V, ("diffML_trig_whtp_d30_g1055_fixedhL40000_L4", 15, tuple([5000]), 4, "whtp", 1.055, 4, "p", 50, 1e-4, None, 2, 4.3))
p.apply_async(PWLD_V, ("pwldML_whtp_var_big_L4", tuple([2000]), 4, "whtp", 1.05, 5, 2.0, 1, "p", 50, 1e-4, None))
p.apply_async(PWLD_V, ("pwldML_whtp_var_small_L4", tuple([2000]), 4, "whtp", 1.05, 5, 1.0/13.0, 1, "p", 50, 1e-4, None))
p.apply_async(PWLD_2D, ("pwld2D_grid2000_var_small_L4", tuple([2000, 2000]), 4, "whtp", 1.06, 5, 1.0/25.0, 1, "p", 50, 1e-3, None))
p.apply_async(PWLD_2D, ("pwld2D_grid2000_var_big_L4", tuple([2000, 2000]), 4, "whtp", 1.06, 5, 2.0, 1, "p", 50, 1e-3, None))



p.close()
p.join()

# Same calls without using the hpc system:
DIFF_V("diffML_trig_whtp_d20_g1055_grid2000", 10, tuple([2000]), 4, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.3)
DIFF_V("diffML_trig_whtp_d30_g1055_grid2000", 15, tuple([2000]), 4, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.3)
DIFF_V("diffML_trig_whtp_d30_g1055_fixedhL40000_L1", 15, tuple([40000]), 1, "whtp", 1.055, 1, "p", 50, 1e-4, None, 2, 4.)
DIFF_V("diffML_trig_whtp_d30_g1055_fixedhL40000_L2", 15, tuple([20000]), 2, "whtp", 1.055, 2, "p", 50, 1e-4, None, 2, 4.))
DIFF_V("diffML_trig_whtp_d30_g1055_fixedhL40000_L3", 15, tuple([10000]), 3, "whtp", 1.055, 3, "p", 50, 1e-4, None, 2, 4.3)
DIFF_V("diffML_trig_whtp_d30_g1055_fixedhL40000_L4", 15, tuple([5000]), 4, "whtp", 1.055, 4, "p", 50, 1e-4, None, 2, 4.3)
PWLD_V("pwldML_whtp_var_big_L4", tuple([2000]), 4, "whtp", 1.05, 5, 2.0, 1, "p", 50, 1e-4, None)
PWLD_V("pwldML_whtp_var_small_L4", tuple([2000]), 4, "whtp", 1.05, 5, 1.0/13.0, 1, "p", 50, 1e-4, None)
PWLD_2D("pwld2D_grid2000_var_small_L3_g107_d25", tuple([2000, 2000]), 3, "whtp", 1.07, 5, 1.0/25.0, 1, "p", 50, 1e-3, None)
PWLD_2D("pwld2D_grid2000_var_big_L3_g107_d25", tuple([2000, 2000]), 3, "whtp", 1.07, 5, 2.0, 1, "p", 50, 1e-3, None)
SPWLD("smallDim6_PWLD_whtp_var_small_grid200_L4_g105", tuple([200]), 4, "whtp", 1.05, 5, 1.0/6.0, 1, "p", 50, 1e-4, 100)
SPWLD("smallDim6_PWLD_whtp_var_big_grid200_L4_g105", tuple([200]), 4, "whtp", 1.05, 5, 2.0, 1, "p", 50, 1e-4, 100)
SPWLD("smallDim6_PWLD_whtp_var_small_grid2000_L4_g105", tuple([2000]), 4, "whtp", 1.05, 5, 1.0/6.0, 1, "p", 50, 1e-4, 100)
SPWLD("smallDim6_PWLD_whtp_var_big_grid2000_L4_g105", tuple([2000]), 4, "whtp", 1.05, 5, 2, 1, "p", 50, 1e-4, 100)
SPWLD("smallDim6_PWLD_whtp_var_big_grid2000_L4_g105_datconstant20", tuple([2000]), 4, "whtp", 1.05, 5, 2, 1, "p", 50, 1e-4, 100, 20)
SPWLD("smallDim6_PWLD_whtp_var_small_grid2000_L3_g105_datconstant20", tuple([2000]), 3, "whtp", 1.05, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 20)
SPWLD("smallDim6_PWLD_whtp_var_small_grid2000_L4_g1055_datconstant20", tuple([2000]), 4, "whtp", 1.055, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 20)
PWLD5("dim5_PWLD_whtp_var_small_grid200_L4_g105", tuple([200]), 4, "whtp", 1.05, 5, 1.0/5.0, 1, "p", 50, 1e-4, 100)
PWLD5("dim5_PWLD_whtp_var_big_grid200_L4_g105", tuple([200]), 4, "whtp", 1.05, 5, 2.0, 1, "p", 50, 1e-4, 100)
PWLD5("dim5_PWLD_whtp_var_small_grid2000_L4_g105", tuple([2000]), 4, "whtp", 1.05, 5, 1.0/5.0, 1, "p", 50, 1e-4, 100)
PWLD5("dim5_PWLD_whtp_var_big_grid2000_L4_g105", tuple([2000]), 4, "whtp", 1.05, 5, 2, 1, "p", 50, 1e-4, 100)
PWLD_2D_small("pwld2D_grid2000_var_small_L3_g107_d16", tuple([2000, 2000]), 3, "whtp", 1.07, 5, 1.0/25.0, 1, "p", 50, 1e-3, None)
PWLD_2D_small("pwld2D_grid2000_var_big_L3_g107_d16", tuple([2000, 2000]), 3, "whtp", 1.07, 5, 2.0, 1, "p", 50, 1e-3, None)
DIFF_P("diffML_trig_whtp_poly_d20_grid2000_L3_datconstant20", 10, tuple([2000]), 3, "whtp", 1.03, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 20)
DIFF_P("diffML_trig_whtp_poly_d30_grid2000_L3_datconstant20", 15, tuple([2000]), 3, "whtp", 1.03, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 20)
DIFF_P("diffML_trig_bp_poly_d40_grid2000_L3_datconstant8", 20, tuple([2000]), 3, "bp", 1.03, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 8)

DIFF_P("diffML_trig_whtp_poly_d30_grid2000_L3_datconstant10_c1015", 15, tuple([2000]), 3, "whtp", 1.015, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 10)

DIFF_P("diffML_trig_whtp_poly_d30_datconstant10_c1015_fixedhL20000_L1", 15, tuple([20000]), 1, "whtp", 1.015, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 10)
DIFF_P("diffML_trig_whtp_poly_d30_datconstant10_c1015_fixedhL20000_L2", 15, tuple([10000]), 2, "whtp", 1.015, 1./2., 2, "p", 50, 1e-3, 100, 2.0, 4.3, 10)
DIFF_P("diffML_trig_whtp_poly_d30_datconstant10_c1015_fixedhL20000_L3", 15, tuple([5000]), 3, "whtp", 1.015, 1./2., 3, "p", 50, 1e-3, 100, 2.0, 4.3, 10)
DIFF_P("diffML_trig_whtp_poly_d30_datconstant10_c1015_fixedhL20000_L4", 15, tuple([2500]), 4, "whtp", 1.015, 1./2., 4, "p", 50, 1e-3, 100, 2.0, 4.3, 10)


DIFF_P("diffML_trig_bp_poly_d40_datconstant8_c102_fixedhL20000_L1", 20, tuple([20000]), 1, "bp", 1.02, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 8)
DIFF_P("diffML_trig_bp_poly_d40_datconstant8_c102_fixedhL20000_L2", 20, tuple([10000]), 2, "bp", 1.02, 1./2., 2, "p", 50, 1e-3, 100, 2.0, 4.3, 8)
DIFF_P("diffML_trig_bp_poly_d40_datconstant8_c102_fixedhL20000_L3", 20, tuple([5000]), 3, "bp", 1.02, 1./2., 3, "p", 50, 1e-3, 100, 2.0, 4.3, 8)
DIFF_P("diffML_trig_bp_poly_d40_datconstant8_c102_fixedhL20000_L4", 20, tuple([2500]), 4, "bp", 1.02, 1./2., 4, "p", 50, 1e-3, 100, 2.0, 4.3, 8)


DIFF_P("diffML_trig_whtp_poly_d40_grid2000_L3_datconstant10_c1015", 20, tuple([2000]), 3, "whtp", 1.015, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 10)
DIFF_P("diffML_trig_bp_poly_d30_grid2000_L4_datconstant10_c1015", 15, tuple([2000]), 4, "bp", 1.015, 1./2., 1, "p", 50, 1e-3, 100, 2.0, 4.3, 10)

SPWLD("smallDim6_PWLD_whtp_var_small_grid2000_L3_g105_datconstant25", tuple([2000]), 3, "whtp", 1.05, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 25)
SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g105_datconstant25", tuple([2000]), 3, "bp", 1.05, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 25)
SPWLD("smallDim6_PWLD_whtp_var_big_grid2000_L3_g105_datconstant25", tuple([2000]), 3, "whtp", 1.05, 5, 2., 1, "p", 50, 1e-4, 100, 25)
SPWLD("smallDim6_PWLD_whtp_var_small_grid2000_L3_g104_datconstant30", tuple([2000]), 3, "whtp", 1.04, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 25)
SPWLD("smallDim6_PWLD_whtp_var_big_grid2000_L3_g104_datconstant30", tuple([2000]), 3, "whtp", 1.04, 5, 2., 1, "p", 50, 1e-4, 100, 25)
SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g105_datconstant25", tuple([2000]), 3, "bp", 1.05, 5, 1.0/6., 1, "p", 50, 1e-4, 100, 25)

SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g1055_datconstant20", tuple([2000]), 3, "bp", 1.055, 5, 1.0/6., 1, "p", 50, 1e-5, 100, 20)


# Some tests for graphing the coefficients nicely
SPWLD("smallDim6_PWLD_womp_var_small_grid2000_L3_g107_datconstant10", tuple([2000]), 3, "womp", 1.07, 5, 1.0/6., 1, "p", 50, 1e-5, 100, 10)
SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g107_datconstant12_HQ_lowTol", tuple([2000]), 3, "bp", 1.07, 5, 1.0/6., 1, "p", 50, 1e-6, 100, 12)
SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g107_datconstant20_HQ_lowTol", tuple([2000]), 3, "bp", 1.07, 5, 1.0/6., 1, "p", 50, 1e-6, 100, 20)
SPWLD("smallDim6_PWLD_bp_var_small_grid2000_L3_g107_datconstant10_HQ_lowTol", tuple([2000]), 3, "bp", 1.07, 5, 1.0/6., 1, "p", 50, 1e-6, 100, 10)

# Make sure the 2D case is easily plotable
PWLD_2D_small("pwld2D_whtp_grid2000_var_small_L3_g107_d16_datconstant8", tuple([2000, 2000]), 3, "whtp", 1.1, 5, 1.0/25.0, 1, "p", 50, 1e-3, 100, 8)
PWLD_2D_small("pwld2D_whtp_grid2000_var_big_L3_g107_d16_datconstant8", tuple([2000, 2000]), 3, "whtp", 1.1, 5, 2.0, 1, "p", 50, 1e-3, 100, 8)
PWLD_2D_small("pwld2D_bp_grid2000_var_small_L3_g107_d16_datconstant8", tuple([2000, 2000]), 3, "bp", 1.1, 5, 1.0/25.0, 1, "p", 50, 1e-3, 100, 8)
PWLD_2D_small("pwld2D_bp_grid2000_var_big_L3_g107_d16_datconstant8", tuple([2000, 2000]), 3, "bp", 1.1, 5, 2.0, 1, "p", 50, 1e-3, 100, 8)
