from dolfin import *

import numpy as np

import shelve

Ntests = 10000

#infile  = 'pwldML_whtp_var_1'
infile  = 'diffML_whtp_d5_g1055p_2000'
# infile  = 'makeSureItWorks'
# Load the results from the MLCSPG 
results = sorted(shelve.open(infile).values(), key=lambda r: r.L)
# for a_result in results:
    # print a_result.L
    # print a_result.wr_model
    # print a_result. cspde_result
first_result = results[0] # At that moment, first_result is a multi-level result
d = first_result.cspde_result[0].d # number of parameters

L_to_check = [r.L for r in results]
print np.max(L_to_check)

Z_tests = np.random.uniform(-1, 1, (Ntests, d))
# Z_tests = first_result.wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (Ntests, d)))


## Make sure we have a good enough accuracy
y_true = first_result.spde_model.samples(Z_tests, 2**(np.max(L_to_check)+2) )
# print y_true

y_recon = []
l2_error = []
l1_error = []
linf_error = []
#l2_error = np.zeros(len(L_to_check))
#l1_error = np.zeros(len(L_to_check))
#linf_error = np.zeros(len(L_to_check))
for cur_results in results:
    cur_recon = cur_results.wr_model.estimate_ML_samples(cur_results.cspde_result, Z_tests)
    l1_error.append(1./Ntests * np.sum(np.abs(y_true - cur_recon)))
    l2_error.append(np.sqrt(1./Ntests * np.sum(np.abs(y_true - cur_recon)**2) ))
    linf_error.append((np.abs(y_true-cur_recon)).max())
    y_recon.append(cur_recon)

for (idx, L) in enumerate(L_to_check):
    print('L = {0} -- L1 error = {1} - L2 error = {2} - Linf error = {3}'.format(L,l1_error[idx],l2_error[idx],linf_error[idx]))