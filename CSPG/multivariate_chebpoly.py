import numpy as np 

# This should compute the tensor product of chebychev polynomial of degrees a_nu (being a multi-index)
def multi_var_chebpoly(a_val, a_nu, with_cheb_weights = True):
    # first check that a_val (the point where to evaluate) and a_nu (the multi-index considered) have the same dimensions
    # Ideally, we would also check that they have the same format and so on... which we are not going to do now! 

    nb_var = len(a_val)
    nb_indices = len(a_nu)
    if nb_var != nb_indices: 
        print('Something went terribly wrong while evaluating multi_var_chebpoly; trying to compute a {0} dimensional polynomial with only {1} coefficients!'.format(nb_var, nb_indices))

    # We should also vectorize this function, i.e. record the number of points at which we wish to estimate the cheb poly
    output_poly = 1
    for idx, value in enumerate(a_val):
        if a_nu[idx] > 0:
            output_poly = output_poly*np.polynomial.chebyshev.chebval(a_val[idx], np.hstack( (np.zeros(a_nu[idx]),1) )) # np.polynomial.chebyshev.chebweight(a_val[idx])
        # elif a_nu[idx] == 1:
            # output_poly = output_poly*np.polynomial.chebyshev.chebval(a_val[idx], 1) # # np.polynomial.chebyshev.chebweight(a_val[idx])

    return output_poly
