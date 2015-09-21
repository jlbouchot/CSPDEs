import numpy as np

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

class Chebyshev(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.sqrt(2)
    name  = 'Chebyshev'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return np.cos(k * np.arccos(x)) * np.sqrt(2)**(k>0)

        return operator_from_matrix(Chebyshev, matrix_from_tensor_indices(J, Z, base, normalization))
