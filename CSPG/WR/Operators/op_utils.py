import numpy as np

#------------------- UTILITIES -------------------#
def chebyshev_table(x, nmax):
    """ Normalised Chebyshev polynomials sqrt(2)^(n>0) T_n(x) for n = 0..nmax, as an array of shape (len(x), nmax+1). """
    T = np.empty((x.size, nmax + 1))
    T[:, 0] = 1.0
    if nmax >= 1:
        T[:, 1] = x
    for n in range(2, nmax + 1):
        T[:, n] = 2.0 * x * T[:, n - 1] - T[:, n - 2]
    T[:, 1:] *= np.sqrt(2.0)
    return T