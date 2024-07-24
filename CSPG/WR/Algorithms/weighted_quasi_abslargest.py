import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2024/07/24"

def weighted_quasi_abslargest(x, s, w):
    '''
    Find the approximate indices of the elements with largest weighted absolute values
        x: input vector
        s: target weighted sparsity
        w: weight sequence
    '''
    sortIndex = np.argsort((np.abs(x)) * (w**(-1)))[::-1]

    k      = 0
    w_test = 0
    w_cur  = 0

    while True:
        w_test = w_cur + w[sortIndex[k]]**2

        if s < w_test:
            break

        w_cur = w_test
        k     = k + 1


    i    = sortIndex[0:k]

    r    = np.zeros_like(x)
    r[i] = x[i]

    return r, i
