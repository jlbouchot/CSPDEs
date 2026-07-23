import numpy as np

import time
import os

from .operator_from_matrix import *
from .LD_bounded_operator  import LD_bounded_operator

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/07/21"

class Chebyshev(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.sqrt(2)
    name  = 'Chebyshev'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return np.cos(k * np.arccos(x)) * np.sqrt(2)**(k>0)

        return operator_from_matrix(Chebyshev, matrix_from_tensor_indices(J, Z, base, normalization))


    def load(data_file, t_fname):
        
        A = operator_from_matrix(Chebyshev, np.load(data_file))
        t_matrix = np.load(t_fname)

        return A, t_matrix

