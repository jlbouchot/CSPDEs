import numpy as np
from scipy.special import lpmn

from operator_from_matrix import *
from LD_bounded_operator  import LD_bounded_operator

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

# TODO: Vectorize for non-scalar z and n
# @np.vectorize # ?
def lpn_normalized(n, z):
    vals, _ = lpmn(0, n, z)

    return vals[0][n] * np.sqrt(n+0.5)


class Legendre(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.inf
    name  = 'Legendre'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return lpn_normalized(k, x) # * 1./np.sqrt(2)

        return operator_from_matrix(Legendre, matrix_from_tensor_indices(J, Z, base, normalization))
