import scipy.integrate as spi
from iterative_solution import compute_true_avg as cta
from iterative_solution import compute_true_avg_alternate as ctaa
from multivariate_chebpoly import multi_var_chebpoly as mvc
from multivariate_chebpoly import vectorized_mvc as vmvc

# from sobol_lib import i4_sobol # for QMC points
# import mcint # for automated MC method

import numpy as np
import math

import os.path


from mpmath import *

def coefs_nquad_wrapper(a_nu): 
    return spi.nquad(lambda x0,x1,x2,x3,x4,x5: ctaa([np.array([x0,x1,x2,x3,x4,x5])], 1.0, 1.0/6.0, 5.0 )*mvc([x0,x1,x2,x3,x4,x5],a_nu),[[-1,1]]*6)

def coefs_mpmath_wrapper(a_nu): # seems like there is no way to use mpmath for more than 3 dimensions :( 
    mp.dps = 10 # 10 digits of accuracy should be more than enough.
    return quad(lambda x0,x1,x2,x3,x4,x5: ctaa([np.array([x0,x1,x2,x3,x4,x5])], 1.0, 1.0/6.0, 5.0 )*mvc([x0,x1,x2,x3,x4,x5],a_nu),[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1], method='tanh-sinh') # tanh-sinh tends to be more accurate when dealing with singularities


# def coefs_dbltpl_wrapper(a_nu): # Will work only for the 6 dimensional case
    # return spi.tplquad(secondIntegrals, -1., 1., fun_mone, fun_one, fun_mone, fun_one, args=a_nu)

def coefs_dbltpl_wrapper(nu0,nu1,nu2,nu3,nu4,nu5): # Will work only for the 6 dimensional case
    return spi.tplquad(secondIntegrals, -1., 1., fun_mone, fun_one, fun_mone, fun_one, args=(nu0,nu1,nu2,nu3,nu4,nu5))

def coefs_MC_wrapper(nb_samples, a_nu, with_cheb_weights = False, rhs = 1., variability = 1./6., abar = 5.): 
    domainsize = 1. # math.pow(2,6)
    np.random.seed(1)
    loc_integrand = lambda x0, x1, x2, x3, x4, x5: integrand(x0,x1,x2,x3,x4,x5, a_nu[0],a_nu[1],a_nu[2],a_nu[3],a_nu[4],a_nu[5], with_cheb_weights, rhs, variability, abar)
    result, error = mcint.integrate(loc_integrand, sampler(), measure=domainsize, n=nb_samples)
    return result

def coefs_QMCSobol_wrapper(nb_samples, a_nu): # note that a lot of calculations will be repeated by doing this. We need to be smarter!
    int_val = 0
    seed = 0
    dim = 6
    for one_sample in xrange(0,nb_samples): 
        [sample, new_seed] = i4_sobol(6, seed)
        int_val = int_val + ctaa(sample, 1.0,1.0/6,5.)*mvc(sample,a_nu)
        seed = new_seed
    return int_val/nb_samples

def faster_QMC_computations(nb_samples, nus): # note that a lot of calculations will be repeated by doing this. We need to be smarter!
    # first find the highest degree considered in the list of nu's
    # then go through the samples, one at a time, and everytime we see a new value, we add to the dictionary together with the T_j(value)
    max_degree = np.max(nus)
    cheb_evals = {}
    weights_eval = {}
    int_val = 0 # as for integral value
    seed = 0 # required for the sobol index
    dim = 6 # that hurts!
    for one_sample in xrange(0,nb_samples): 
        [sample, new_seed] = i4_sobol(6, seed)
        sample = sample*2 -1 # To set the anchor at (-1,-1,-1 ...) instead of the usual (0,0,...) for QMC methods
        not_computed = [a_param for a_param in one_sample if a_param not in cheb_evals] # contains the values that we have not precomputed before
        for to_compute in not_computed:
            # add these guys to the dictionary. 
            to_add_to_dict = [1]
            for one_deg in xrange(1,max_degree+1): 
                to_add_to_dict.append(np.polynomial.chebyshev.chebval(to_compute, np.hstack( (np.zeros(one_deg),1) )))
            cheb_evals[to_compute] = to_add_to_dict
            weights_eval[to_compute] = np.polynomial.chebyshev.chebweight(to_compute)
        int_val = int_val + ctaa(sample-1, 1.0,1.0/6,5.)*mvc(sample,a_nu)
        seed = new_seed
    return int_val/nb_samples

def coefs_Smolyak_wrapper():
    return 1

def integrand(x0,x1,x2,x3,x4,x5,nu0,nu1,nu2,nu3,nu4,nu5, with_weights = True, rhs = 1., variability = 1.0/6., abar = 5.):
    # print ctaa([np.array([x0,x1,x2,x3,x4,x5])], rhs, variability, abar )
    return ctaa([np.array([x0,x1,x2,x3,x4,x5])], rhs, variability, abar )*mvc([x0,x1,x2,x3,x4,x5],np.array([nu0,nu1,nu2,nu3,nu4,nu5]), with_weights)
# def integrand(x0,x1,x2,x3,x4,x5,a_nu):
    # return ctaa([np.array([x0,x1,x2,x3,x4,x5])], 1.0, 1.0/6.0, 5.0 )*mvc([x0,x1,x2,x3,x4,x5],a_nu)

# def secondIntegrals(x3,x4,x5,a_nu):
    # res, err = spi.tplquad(integrand, -1., 1., fun_mone, fun_one, fun_mone, fun_one, args=(x3,x4,x5,a_nu))
    # return res

def secondIntegrals(x3,x4,x5,nu0,nu1,nu2,nu3,nu4,nu5):
    res, err = spi.tplquad(integrand, -1., 1., fun_mone, fun_one, fun_mone, fun_one, args=(x3,x4,x5,nu0,nu1,nu2,nu3,nu4,nu5))
    return res

def mpsecondIntegrals(x3,x4,x5,nu0,nu1,nu2,nu3,nu4,nu5):
    res, err = tplquad(integrand, [-1., 1.], [-1., 1.], [-1., 1.], args=(x3,x4,x5,nu0,nu1,nu2,nu3,nu4,nu5))
    return res

def fun_one(p1 = 0, p2 = 0): # the boundaries of the tplquad have to be functions of 0, 1, or 2 variables!
    return 1.0

def fun_mone(p1 = 0, p2 = 0): # the boundaries of the tplquad have to be functions of 0, 1, or 2 variables!
    return -1.0

def sampler(): # for the Monte-Carlo approach
    while True: # Generates the chebyshev measures via a uniform measure
        x0     = math.cos(np.random.uniform(0,math.pi))
        x1     = math.cos(np.random.uniform(0,math.pi))
        x2     = math.cos(np.random.uniform(0,math.pi))
        x3     = math.cos(np.random.uniform(0,math.pi))
        x4     = math.cos(np.random.uniform(0,math.pi))
        x5     = math.cos(np.random.uniform(0,math.pi))
        yield (x0, x1, x2, x3, x4, x5)

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
            ys = ctaa(zs, rhs, variability, mean_field)
            cheb_vals = vmvc(zs,a_nu)
            final_avg = final_avg + sum(ys*cheb_vals)
            dict_results[hashable_nu][nb_tests] = final_avg/nb_tests
            np.save(dump_fname, dict_results) 

    # one last batch of tests to make sure we have the required number
    remaining_tests = tot_tests - nb_tests
    if remaining_tests > 0:
        nb_tests = nb_tests+remaining_tests
        zs = np.cos(np.random.uniform(0,np.pi,[remaining_tests, 6]))
        ys = ctaa(zs, rhs, variability, mean_field)
        cheb_vals = vmvc(zs,a_nu)
        final_avg = final_avg + sum(ys*cheb_vals)
        dict_results[hashable_nu][nb_tests] = final_avg/nb_tests
        np.save(dump_fname, dict_results) 

    return final_avg/nb_tests
