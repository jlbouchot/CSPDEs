from test_diff_avg_v_ML import Main as DIFF_V
from test_piecewiselineardiffusion_v import Main as PWLD_V # As in PieceWiseLinearDiffusion
from test_piecewiselineardiffusion_v_2D import Main as PWLD_2D
from test_diff_avg_v_ML_2D import Main as DIFF_V_2D

from multiprocessing import Pool


# Compute solutions
p = Pool(4)

p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d20_g1055_grid2000", alpha = 2, d = 10, grid_points = 2000, L_max = 4, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d30_g1055_grid2000", alpha = 2, d = 15, grid_points = 2000, L_max = 4, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d30_g1055_fixedhL40000_L1", alpha = 2, d = 15, grid_points = 40000, L_min = 1, L_max = 1, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d30_g1055_fixedhL40000_L2", alpha = 2, d = 15, grid_points = 20000, L_min = 2, L_max = 2, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d30_g1055_fixedhL40000_L3", alpha = 2, d = 15, grid_points = 10000, L_min = 3, L_max = 3, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(DIFF_V, (outfile = "diffML_trig_whtp_d30_g1055_fixedhL40000_L4", alpha = 2, d = 15, grid_points = 5000, L_min = 4, L_max = 4, algo_name = "whtp", gamma = 1.055, nb_iter = 50, epsilon = 1e-4))
p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_big_L4", grid_points = 2000, L_max = 4, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 2.0, nb_iter = 50, epsilon = 1e-4)
p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_small_L4", grid_points = 2000, L_max = 4, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 1.0/13.0, nb_iter = 50, epsilon = 1e-4)
p.apply_async(PWLD_2D, (outfile = "pwld2D_grid2000_var_small_L4", L_max = 4, algo_name = "whtp", grid_points = tuple([2000, 2000]), gamma = 1.06, abar = 5, variability = 1.0/25.0, nb_iter = 50, epsilon = 1e-3))
p.apply_async(PWLD_2D, (outfile = "pwld2D_grid2000_var_big_L4", L_max = 4, algo_name = "whtp", grid_points = tuple([2000, 2000]), gamma = 1.06, abar = 5, variability = 2.0, nb_iter = 50, epsilon = 1e-3))



p.close()
p.join()