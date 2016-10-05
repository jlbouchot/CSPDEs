# This is a little function that calls FEniCS given the tests made previously.
# I was sick of having to look at what I had done before!

import shelve
import numpy as np

def compute_numerical_solutions(infile, nb_tests, refine_ratio = None):
    print("Loading {0} ...".format(infile))
    results = sorted(shelve.open(infile).values(), key=lambda r: r.L)

    finest_result = results[-1] # At that moment, finest_result is a multi-level result containing the finest information
    d             = finest_result.cspde_result[0].d # number of parameters
    spde_model    = finest_result.spde_model
    epsilon       = finest_result.epsilon 
    # First build the global set of active chebychef polynomials
    nb_lvl = finest_result.L # Corresponds to the number of levels for the finest model

    Z = finest_result.wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (nb_tests, d)))

    if refine_ratio is not None: 
        finest_result.spde_model.refine_mesh(refine_ratio)
    y_new = finest_result.spde_model.samples(Z)
    return y_new, Z, finest_result
