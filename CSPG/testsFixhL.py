from test_diff_avg_v_ML import Main as DIFF_V

from multiprocessing import Pool


# Compute solutions
p = Pool(5)
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L5_fixhL", 5, 2000, 5, "whtp", 1.055, 5))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L4_fixhL", 5, 4000, 5, "whtp", 1.055, 4))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L3_fixhL", 5, 8000, 5, "whtp", 1.055, 3))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L2_fixhL", 5, 16000, 5, "whtp", 1.055, 2))
p.apply_async(DIFF_V, ("diffML_whtp_d5_g1055_L1_fixhL", 5, 32000, 5, "whtp", 1.055, 1))

p.close()
p.join()
