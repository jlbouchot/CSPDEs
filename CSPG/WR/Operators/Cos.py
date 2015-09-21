import numpy as np

from operator_from_matrix import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

class Cos:
    # L_inf norm of the basis functions associated to this operator
    theta = 1.0
    name  = 'Cos'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return (1./np.sqrt(2))**np.count_nonzero(k==0) * np.cos(np.pi * x * k)

        return operator_from_matrix(Cos, matrix_from_tensor_indices(J, Z, base, normalization))
