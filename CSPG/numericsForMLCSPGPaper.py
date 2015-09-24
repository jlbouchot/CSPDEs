from test_diff_avg_v_ML import Main as DIFF_V
from test_piecewiselineardiffusion_v import Main as PWLD_V # As in PieceWiseLinearDiffusion

## Still have to deal with the plot scripts
# from plot_comparison_1D import plot_1D
# from plot_comparison_2D import plot_2D
# from plot_s             import plot_s
# from plot_qmc           import plot_qmc

from multiprocessing import Pool


# Compute solutions
p = Pool(4)
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d20_g105", d = 20, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g105", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g104", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.04))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g103", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.03))
# p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_1", grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 1.0/13.0))
# p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_5", grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 5.0/13.0))
p.apply_async(DIFF_V, ("diffML_whtp_d20_g105", 20, 2000, 5, "whtp", 1.05))
p.apply_async(DIFF_V, ("diffML_whtp_d10_g105", 10, 2000, 5, "whtp", 1.05))
p.apply_async(DIFF_V, ("diffML_whtp_d10_g104", 10, 2000, 5, "whtp", 1.04))
p.apply_async(DIFF_V, ("diffML_whtp_d10_g103", 10, 2000, 5, "whtp", 1.03))
# p.apply_async(PWLD_V, ("pwldML_whtp_var_1", 2000, 5, "whtp", 1.05, 5, 1.0/13.0))
# p.apply_async(PWLD_V, ("pwldML_whtp_var_5", 2000, 5, "whtp", 1.05, 5, 5.0/13.0))

p.close()
p.join()

# # Plot
# plot_1D("r_gamma_diff", "diff_1d", [35, 65, 140], [1.06], [0, 1, 2, 3])
# plot_1D("r_gamma_fin" , "fin_1d" , [35, 95, 155], [1.08], [0, 1, 2, 3])
# plot_1D("r_gamma_cd"  , "cd_1d"  , [35, 95, 155], [1.01+0.005+0.005], [0, 1])

# plot_2D("r_gamma_diff", "diff_2d_01",  95,  1.06, 0, 1)
# plot_2D("r_gamma_diff", "diff_2d_02",  95,  1.06, 0, 2)
# plot_2D("r_gamma_fin" , "fin_2d_01" ,  95,  1.08, 0, 1)
# plot_2D("r_gamma_fin" , "fin_2d_12" ,  95,  1.08, 1, 2)
# plot_2D("r_gamma_cd"  , "cd_2d_s35"     ,  35, 1.01+0.005+0.005, 0, 1)
# plot_2D("r_gamma_cd"  , "cd_2d_s155"     ,  155, 1.01, 0, 1)

# plot_s("r_gamma_diff", "diff_s")
# plot_s("r_gamma_fin" , "fin_s")
# plot_s("r_gamma_cd"  , "cd_s")

# plot_qmc("r_gamma_diff", "diff")
# plot_qmc("r_gamma_fin" , "fin")
# plot_qmc("r_gamma_cd"  , "cd")

