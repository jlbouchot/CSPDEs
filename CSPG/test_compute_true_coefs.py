from compare_true_coefs import get_true_coefs as gtc
from compare_true_coefs import get_computed_coefs as gcc

import numpy as np

true_coefs, J = gtc("smallTestsPoisson", None)
print type(J)
estimations, other_J = gcc("smallTestsPoisson", None, J) # other_J should be the same as J itself
print type(other_J)
print type(true_coefs)
print type(estimations)
sorted_true = np.sort(np.abs(true_coefs))[::-1]
sorted_estimations = np.sort(np.abs(estimations.tolist()))[::-1]
print sorted_true[0:10]
print sorted_estimations[0:10]

J = np.asarray(J)
other_J = np.asarray(other_J)

print type(J)
print type(other_J)

print J.shape
print other_J.shape

d = 10
a_nu = [0,0,0,0,0,0,0,0,0,0]
another_nu = [1,0,0,0,0,0,0,0,0,0]
find_idx = np.array(np.sum(J == a_nu, axis = 1))
cur_idx_true = find_idx.tolist().index(d)
find_idx = np.array(np.sum(other_J == a_nu, axis = 1))
cur_idx_estimations = find_idx.tolist().index(d)
print("First coefficient estimated as position {0} and with value {1} in the true coefficients, while it has position {2} and value {3} in the estimated coefficients".format(cur_idx_true, sorted_true[cur_idx_true], cur_idx_estimations, sorted_estimations[cur_idx_estimations]))

find_idx = np.array(np.sum(J == another_nu, axis = 1))
cur_idx_true = find_idx.tolist().index(d)
find_idx = np.array(np.sum(other_J == another_nu, axis = 1))
cur_idx_estimations = find_idx.tolist().index(d)
print("Second coefficient estimated as position {0} and with value {1} in the true coefficients, while it has position {2} and value {3} in the estimated coefficients".format(cur_idx_true, sorted_true[cur_idx_true], cur_idx_estimations, sorted_estimations[cur_idx_estimations]))

print("Second coefficient estimated as position {0} and with value {1} in the true coefficients, while it has position {2} and value {3} in the estimated coefficients".format(cur_idx_true, sorted_true[cur_idx_true]/sorted_true[0]*sorted_estimations[0], cur_idx_estimations, sorted_estimations[cur_idx_estimations]))

print np.linalg.norm(sorted_estimations-sorted_true/sorted_true[0]*sorted_estimations[0])
print np.linalg.norm(sorted_estimations-sorted_true)

# blah = tuple([np.array([1,0,0,0]), np.array([0,0,0,0]), np.array([0,1,0,0])])
