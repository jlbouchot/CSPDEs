from test_diff_avg_v_ML_2D import Main as DIFF_V_2D

## Still have to deal with the plot scripts
# from plot_comparison_1D import plot_1D
# from plot_comparison_2D import plot_2D
# from plot_s             import plot_s
# from plot_qmc           import plot_qmc


# Compute solutions
# p = Pool(4)
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d20_g105", d = 20, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g105", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g104", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.04))
# p.apply_async(DIFF_V, (outfile = "diffML_whtp_d10_g103", d = 10, grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.03))
# p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_1", grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 1.0/13.0))
# p.apply_async(PWLD_V, (outfile = "pwldML_whtp_var_5", grid_points = 2000, L_max = 5, algo_name = "whtp", gamma = 1.05, abar = 5, variability = 5.0/13.0))


DIFF_V_2D(outfile = "first2Dtests", d = 5, L_max = 2, algo_name = "whtp", grid_points = tuple([1000, 1000]), gamma = 1.06, nb_iter = 50, nb_tests = 50)
