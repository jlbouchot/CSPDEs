import numpy as np

from operator_from_matrix import *
from LD_bounded_operator  import LD_bounded_operator


class Chebyshev(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.sqrt(2)
    name  = 'Chebyshev'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return np.cos(k * np.arccos(x)) * np.sqrt(2)**(k>0)

        return operator_from_matrix(Chebyshev, matrix_from_tensor_indices(J, Z, base, normalization))
