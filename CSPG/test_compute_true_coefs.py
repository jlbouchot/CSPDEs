from compare_true_coefs import get_true_coefs as gtc
from compare_true_coefs import get_computed_coefs as gcc

import numpy as np

true_coefs, J = gtc("smallTestsPoisson", None)
print true_coefs
print len(true_coefs)
print len(J)
estimations = gcc("smallTestsPoisson", None, J)
print type(true_coefs)
print type(estimations)
sorted_true = np.sort(np.abs(true_coefs))[::-1]
sorted_estimations = np.sort(np.abs(estimations.tolist()))[::-1]
print sorted_true[0:10]
print sorted_true[0:10]
