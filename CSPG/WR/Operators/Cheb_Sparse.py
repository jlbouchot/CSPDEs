"""
Matrix-free tensorised Chebyshev sampling operator exploiting sparse multi-indices.

Same maths, same interface and same results (to rounding) as Cheb_Alt, but faster:

1. Sparse multi-indices. Entry (i, j) of the sampling matrix is
       A[i, j] = (1/norm) * prod_{k=1}^{d} T_{nu_j,k}(z_{i,k})      (normalised T_n, T_0 = 1)
   Cheb_Alt multiplies all d factors, although nu_j has only a handful of non-zero
   entries (|supp nu| <= ~3 for the weighted index sets J_s). Here each nu_j is stored
   as the list of its *non-zero* entries only, so the cost of a product is
   O(m * sum_j |supp nu_j|) instead of O(m * N * d): typically 5-10x less work for d = 20.

2. Sparse iterates. In WIHT/WHTP the vector x is s-sparse. apply() only loops over
   the non-zero entries of x: O(m * |supp x|) instead of O(m * N).

3. Real multi-threading. operator_from_matrix_Alt uses prange inside @njit *without*
   parallel=True, so it runs on one core. Here the outer loops are parallel
   (numba threads, set NUMBA_NUM_THREADS to control).

4. Precomputed univariate table. The (m x C) table of T_n(z_{i,k}) is built with a
   vectorised recurrence instead of a triple Python loop calling cos(k arccos x).

5. to_dense(): when m*N*8 bytes fits in memory, materialising A once and using BLAS
   (multi-threaded dgemv) beats any matrix-free scheme over hundreds of iterations.

Usage: replace WR.Operators.Cheb_Alt by WR.Operators.Cheb_Sparse in the driver.
"""

import os
import numpy as np
from numba import njit, prange

from .LD_bounded_operator import LD_bounded_operator
from .op_utils import chebyshev_table
from .operator_from_sparse_tensor import operator_from_sparse_tensor, sparse_tensor_from_indices, _npz

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2026, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/10/14"
__created__ = "2026/10/09"

class Cheb_Sparse(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.sqrt(2)
    name = 'Cheb_Sparse'

    @staticmethod
    def create(J, Z, normalization=None):
        T, ptr, cols, scale = sparse_tensor_from_indices(J, Z, chebyshev_table, normalization)
        return operator_from_sparse_tensor(Cheb_Sparse, T, ptr, cols, scale, J)

    @staticmethod
    def load(data_file, t_fname):
        f = np.load(_npz(data_file))
        return operator_from_sparse_tensor(Cheb_Sparse, f["T"], f["ptr"], f["cols"], f["scale"]), np.load(t_fname)
