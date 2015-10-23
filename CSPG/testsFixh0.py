from test_diff_avg_v_ML import Main as DIFF_V

from multiprocessing import Pool


# Compute solutions
p = Pool(5)
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L5", 5, 2000, 5, "whtp", 1.055, 5))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L4", 5, 2000, 5, "whtp", 1.055, 4))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L3", 5, 2000, 5, "whtp", 1.055, 3))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L2", 5, 2000, 5, "whtp", 1.055, 2))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L1", 5, 2000, 5, "whtp", 1.055, 1))

p.close()
p.join()
