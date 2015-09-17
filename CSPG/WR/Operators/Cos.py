import numpy as np

from operator_from_matrix import *

class Cos:
    # L_inf norm of the basis functions associated to this operator
    theta = 1.0
    name  = 'Cos'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return (1./np.sqrt(2))**np.count_nonzero(k==0) * np.cos(np.pi * x * k)

        return operator_from_matrix(Cos, matrix_from_tensor_indices(J, Z, base, normalization))
