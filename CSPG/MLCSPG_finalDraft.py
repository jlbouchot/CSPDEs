from test_diff_avg_v_ML import Main as DIFF_V
from test_piecewiselineardiffusion_v import Main as PWLD_V # As in PieceWiseLinearDiffusion
from test_piecewiselineardiffusion_v_2D import Main as PWLD_2D
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
