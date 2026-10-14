"""
Cheb_Alt++: same operator as Cheb_Alt, but the univariate tensors are built with the
three-term recurrence
    T_0(x) = 1,  T_1(x) = x,  T_{n+1}(x) = 2 x T_n(x) - T_{n-1}(x)
vectorised over the samples, instead of a triple Python loop calling cos(n arccos x).
The resulting operator (operator_from_matrix_Alt) and the save/load format are those of Cheb_Alt.
"""

import numpy as np
import os.path

from .operator_from_matrix_Alt import operator_from_matrix_Alt, matrix_from_tensor_indices, univ_tensor_from_table_indices
from .op_utils import chebyshev_table
from .LD_bounded_operator import LD_bounded_operator

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2026, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/09"
__lastmodified__ = "2026/10/09"

class Cheb_AltPP(LD_bounded_operator):
    theta = np.sqrt(2)
    name  = 'Cheb_AltPP'

    @staticmethod
    def create(J, Z, normalization=None):
        return operator_from_matrix_Alt(Cheb_AltPP, matrix_from_tensor_indices(J, Z, None, normalization), univ_tensor_from_table_indices(J, Z, chebyshev_table, normalization), J)

    @staticmethod
    def load(data_file, t_fname):
        dirname = os.path.dirname(data_file)
        bname = os.path.basename(data_file)
        A = np.load(os.path.join(dirname, "A_"+bname))
        # Univariate blocks are saved as "{idx}_univariate_<bname>.npy", sorted by their integer prefix
        univariate_files = [f for f in os.listdir(dirname) if f.endswith("_univariate_"+bname)]
        univariate_files.sort(key=lambda s: int(s.split("_univariate_")[0]))
        univariate = [np.load(os.path.join(dirname, fname)) for fname in univariate_files]
        J = np.load(os.path.join(dirname, "J_"+bname))

        return operator_from_matrix_Alt(Cheb_AltPP, A, univariate, J), np.load(t_fname)
