import numpy as np

from dolfin import *

import sys
import shelve


def get_true_coefs(infile, outfile, J = None):
    # We basically oversample over the parameter space and compute the L2 minimization of the coefficients in the union of all J_s
    print("Loading {0} ...".format(infile))
    results = sorted(shelve.open(infile).values(), key=lambda r: r.L)

    ### Plot
    finest_result = results[-1] # At that moment, finest_result is a multi-level result containing the finest information
    d             = finest_result.cspde_result[0].d # number of parameters
    spde_model    = finest_result.spde_model
    epsilon       = finest_result.epsilon 
    # First build the global set of active chebychef polynomials
    nb_lvl = finest_result.L # Corresponds to the number of levels for the finest model

    # Build the seet of potential active components, i.e. J = cup for all level J_l (that is, if not already passed as an input)
    if J is None:
        J = []
        J.append(np.zeros(d, dtype=int))

        for one_lvl in range(nb_lvl):
            J_s = finest_result.cspde_result[one_lvl].J_s
            for a_nu in J_s:
                if not np.any(np.sum(J == a_nu, axis = 1) == d):
                    J.append(a_nu)
    

    # Now that we have the full 'interesting set' we can oversample it and use an L2 minimization
    N = len(J)
    m = int(1.1*N) # Why 10? Well... why not?
    # Compute the right hand side: 
    print("Computing {0} samples to estimate {1} coefficients".format(m,N))
    Z = finest_result.wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
    y_recon = finest_result.wr_model.estimate_ML_samples(finest_result.cspde_result, Z)
    y_new = finest_result.spde_model.samples(Z)
    # Compare
    difference = np.abs(y_new - y_recon)
    print("Maximum error: {0}; Average error: {1} taken from {2} sample points".format(difference.max(), difference.sum()/m, m))
    # Compute the 'sensing matrix' to be pseudo inverted
    print("Creating matrix before inversion -- this *will* take some time!")
    A = finest_result.wr_model.operator.create(J, Z, normalization=1).A
    # A = finest_result.wr_model.operator.create(J, Z, normalization=np.sqrt(m)).A
    print("Solving linear system -- this may take some time")
    # recon_coefs = np.linalg.lstsq(A,y_recon)
    true_coefs = np.linalg.lstsq(A,y_new)
    # print("Solving linear system with pinv -- this may take some time")
    # pinv_coefs = np.linalg.lstsq(A,y_new, rcond=1e-10)
    print("Norm of the error between true values and l2 optimized:{0}".format(np.linalg.norm(y_new-np.dot(A,true_coefs[0]))))
    # print("Norm of the error between true values and l2 optimized (with reconstructed y's):{0}".format(np.linalg.norm(y_recon-np.dot(A,true_coefs[0]))))
    # print("Norm of the error between true values and pinv results:{0}".format(np.linalg.norm(y_new-np.dot(A,pinv_coefs[0]))))




    # Still have to write our findings into an output file

    return true_coefs[0], J


def get_computed_coefs(infile,outfile, J = None):
    # Have to go through all the single levels (at the finest goal accuracy) and add them alltogether
    print("Computing the final multi-level estimation of the coefficients")
    print("Loading {0} ...".format(infile))
    results = sorted(shelve.open(infile).values(), key=lambda r: r.L)

    ### Plot
    finest_result = results[-1] # At that moment, finest_result is a multi-level result containing the finest information
    d             = finest_result.cspde_result[0].d # number of parameters
    spde_model    = finest_result.spde_model
    epsilon       = finest_result.epsilon 
    # First build the global set of active chebychef polynomials
    nb_lvl = finest_result.L # Corresponds to the number of levels for the finest model

    # First, build the global active set (as in the function above)    
    if J is None:
        J = []
        J.append(np.zeros(d, dtype=int))

        for one_lvl in range(nb_lvl):
            J_s = finest_result.cspde_result[one_lvl].J_s
            for a_nu in J_s:
                if not np.any(np.sum(J == a_nu, axis = 1) == d):
                    J.append(a_nu)

    N = len(J)

    estimated_coefs = np.zeros(N)
    for idx_nu, a_nu in enumerate(J): 
        # Check all the coefficients one after the other, and see in which level they appear
        for one_lvl in range(nb_lvl): 
            find_idx = np.array(np.sum(finest_result.cspde_result[one_lvl].J_s == a_nu, axis = 1))
            if np.any(find_idx == d):
                cur_idx = find_idx.tolist().index(d)
                # print type(estimated_coefs[idx_nu])
                # print type(finest_result.cspde_result[one_lvl].result.x[cur_idx])
                # print("Index for the coefficients {0} : {1}. It has magnitude {2} at level {3}".format(a_nu, cur_idx, finest_result.cspde_result[one_lvl].result.x[cur_idx], one_lvl))
                estimated_coefs[idx_nu] = estimated_coefs[idx_nu] + finest_result.cspde_result[one_lvl].result.x[cur_idx]/np.sqrt(finest_result.cspde_result[one_lvl].m) 

    return estimated_coefs, J
