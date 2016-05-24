from test_piecewiselineardiffusion_v_2D import Main as PWLD_2D


PWLD_2D(outfile = "PWLD2D_firstTests", L_max = 2, algo_name = "whtp", grid_points = tuple([200, 200]), gamma = 1.06, abar = 5, variability = None, nb_iter = 50, nb_tests = 50, epsilon = 1e-3)
