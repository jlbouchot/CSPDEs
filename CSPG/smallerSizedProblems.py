from test_diff_avg_v_ML import Main as DIFF_V
from test_piecewiselineardiffusion_v import Main as PWLD_V # As in PieceWiseLinearDiffusion

from multiprocessing import Pool


# Compute solutions
p = Pool(3)
p.apply_async(DIFF_V, ("diffML_whtp_d6_g105", 6, 2000, 5, "whtp", 1.05))
p.apply_async(DIFF_V, ("diffML_whtp_d6_g105", 6, 2000, 5, "whtp", 1.04))
p.apply_async(DIFF_V, ("diffML_whtp_d6_g104", 6, 2000, 5, "whtp", 1.03))

p.close()
p.join()
