import numpy as np
import os
from multivariate_chebpoly import vectorized_mvc as vmvc



def compute_true_avg(ys, rhs, variability, mean_field):
    # print type(ys)
    # print ys.shape
    if type(ys) is list:
        ys = np.array(ys) 
    if len(ys.shape) == 2:
        nb_ys, nb_param = ys.shape
    else: 
        nb_param = ys.shape[0]
        nb_ys = 1

    delta_ts = 1.0/(nb_param) # length of the segments

    one_over_a = np.divide(1.0,mean_field+variability*ys)
    harm_means = np.cumsum(one_over_a,axis = 1)
    # weighted_means = np.cumsum(np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys).transpose()*one_over_a,axis = 1)
    weighted_means = np.cumsum(np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys)*one_over_a,axis = 1)

    if nb_ys == 1: 
        C_2 = np.divide(1.0,harm_means[-1])/delta_ts + rhs*delta_ts / 2.0 * np.divide(weighted_means[-1],harm_means[-1])
        C_2 = C_2[-1] # This is really really dirty! We should improve on this!
        # u_tilda_k = delta_ts*C_2*one_over_a - rhs*delta_ts**2/2*np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys)*one_over_a
        u_tilda_k = delta_ts*C_2*one_over_a - rhs*delta_ts**2/2*np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys)*one_over_a
        u_tilda_k = np.hstack(([[0.0]], u_tilda_k ) )
        u_k = np.cumsum(u_tilda_k, axis = 1)
        to_add = delta_ts**2/2*C_2*one_over_a - rhs/6.0*delta_ts**3*one_over_a*np.array([3.0*np.linspace(1,nb_param,nb_param)-2.0]*nb_ys)
        return_val = np.sum(delta_ts*u_k[:,0:-1] + to_add, axis = 1)
    else: 
        C_2 = np.divide(1.0,harm_means[:,-1])/delta_ts + rhs*delta_ts / 2.0 * np.divide(weighted_means[:,-1],harm_means[:,-1])
        u_tilda_k = delta_ts*np.array([C_2]*nb_param).transpose()*one_over_a - rhs*delta_ts**2/2*np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys)*one_over_a
        u_tilda_k = np.hstack((np.zeros((nb_ys,1), dtype=float), u_tilda_k ) )
        u_k = np.cumsum(u_tilda_k, axis = 1)
        to_add = delta_ts**2/2*np.array([C_2]*nb_param).transpose()*one_over_a - rhs/6.0*delta_ts**3*one_over_a*np.array([3.0*np.linspace(1,nb_param,nb_param)-2.0]*nb_ys)
        return_val = np.sum(delta_ts*u_k[:,0:-1] + to_add, axis = 1)
    

    return return_val




def vectorized_MC_integration(tot_tests, a_nu, batch_size = 100000, rhs = 1., variability = 1./6., mean_field = 5., regular_saves = True, dump_fname = None): # This can be adapted to make it better for later tests. At the moment, it is implemented fast and dirty for our particular setup

    final_avg = 0

    # Since lists, nd_arrays are not hashable, we need to create a key (unique) for each nu... 
    hashable_nu = '{}'.format(a_nu)
    print('******************* Current nu = ' + hashable_nu)

    ## This definitely needs to be done better! 
    if dump_fname is None:
        dump_fname = 'MCEstimations_PWLD_6dim_small_var.npy'
    if os.path.isfile(dump_fname) is False:
        dict_results = {hashable_nu: {}} # or load if the file already exists
        batches = range(batch_size, tot_tests + 1, batch_size)
    else:
        dict_results = np.load(dump_fname).item()
        if hashable_nu not in dict_results:
            dict_results[hashable_nu] = {}
            nb_tests = 0
            batches = range(batch_size, tot_tests + 1, batch_size)
        else: 
            computed_tests = sorted(dict_results[hashable_nu].keys())
            if tot_tests >= computed_tests[-1]:
                batches = range(computed_tests[-1]+batch_size, tot_tests+1, batch_size)
                nb_tests = computed_tests[-1]
                final_avg = dict_results[hashable_nu][computed_tests[-1]]*computed_tests[-1]
            else:
                nb_tests = computed_tests[[computed_tests[i]-tot_tests >= 0 for i in xrange(0,len(computed_tests))].index(True)-1]
                batches = range(nb_tests+batch_size, tot_tests+1, batch_size)
                final_avg = dict_results[hashable_nu][nb_tests]*nb_tests


    for a_batch in batches:
        print('\tCurrent number of test samples: {}'.format(a_batch)) 
        if a_batch in dict_results[hashable_nu]: 
            # Don't redo the calculation, and load the one already done:
            nb_tests = a_batch
            final_avg = dict_results[hashable_nu][nb_tests]*nb_tests
        else: 
            nb_tests = a_batch
            zs = np.cos(np.random.uniform(0,np.pi,[batch_size, 6]))
            ys = compute_true_avg(zs, rhs, variability, mean_field)
            cheb_vals = vmvc(zs,a_nu)
            final_avg = final_avg + sum(ys*cheb_vals)
            dict_results[hashable_nu][nb_tests] = final_avg/nb_tests
            np.save(dump_fname, dict_results) 

    # one last batch of tests to make sure we have the required number
    remaining_tests = tot_tests - nb_tests
    if remaining_tests > 0:
        nb_tests = nb_tests+remaining_tests
        zs = np.cos(np.random.uniform(0,np.pi,[remaining_tests, 6]))
        ys = compute_true_avg(zs, rhs, variability, mean_field)
        cheb_vals = vmvc(zs,a_nu)
        final_avg = final_avg + sum(ys*cheb_vals)
        dict_results[hashable_nu][nb_tests] = final_avg/nb_tests
        np.save(dump_fname, dict_results) 

    return final_avg/nb_tests
