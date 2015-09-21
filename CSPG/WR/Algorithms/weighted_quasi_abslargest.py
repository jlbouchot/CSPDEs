import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

def weighted_quasi_abslargest(x, s, w):
    # Get indices of elements sorted in descending order
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
