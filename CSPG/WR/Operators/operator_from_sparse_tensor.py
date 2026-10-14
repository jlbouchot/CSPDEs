import numpy as np
import os.path
from numba import njit,prange

from .op_utils import chebyshev_table

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/10/14"
__created__ = "2026/10/09"

class operator_from_sparse_tensor:
    '''
    CSR representation of operators
    '''
    def __init__(self, model, T, ptr, cols, scale, J=None):
        self.model = model
        self.T = np.ascontiguousarray(T)          # (m, C): univariate values, one column per (dim, degree)
        self.TT = np.ascontiguousarray(T.T)       # (C, m): same, contiguous in the sample index for the adjoint
        self.ptr = np.asarray(ptr, dtype=np.int64)
        self.cols = np.asarray(cols, dtype=np.int64)
        self.scale = float(scale)
        self.m = self.T.shape[0]
        self.n = self.ptr.size - 1
        self.J = J
        self._dense = None

    @property
    def A(self):
        return self.to_dense()

    def to_dense(self, max_bytes=None):
        if self._dense is None:
            if max_bytes is not None and 8 * self.m * self.n > max_bytes:
                raise MemoryError("Dense sampling matrix would need {0:.1f} GB".format(8e-9 * self.m * self.n))
            self._dense = _columns(self.T, self.ptr, self.cols, np.arange(self.n, dtype=np.int64), self.scale)
        return self._dense

    def apply(self, x):
        if self._dense is not None:
            return self._dense @ x
        idx = np.flatnonzero(x).astype(np.int64)
        return apply_csr(self.T, self.ptr, self.cols, idx, np.ascontiguousarray(x[idx], dtype=np.float64), self.scale)

    def apply_adj(self, y):
        if self._dense is not None:
            return self._dense.T @ y
        return apply_adj_csr(self.TT, self.ptr, self.cols, np.ascontiguousarray(y, dtype=np.float64), self.scale)

    def genSubMatrix(self, U):
        U = np.asarray(U, dtype=np.int64)
        if self._dense is not None:
            return self._dense[:, U]
        return _columns(self.T, self.ptr, self.cols, U, self.scale)

    def save(self, data_file):
        np.savez(_npz(data_file), T=self.T, ptr=self.ptr, cols=self.cols, scale=self.scale)

@njit(parallel=True, cache=True)
def apply_csr(T, ptr, cols, x_idx, x_val, scale):
    m = T.shape[0]
    res = np.empty(m)
    for i in prange(m):
        s = 0.0
        for jj in range(x_idx.size):
            j = x_idx[jj]
            p = x_val[jj]
            for q in range(ptr[j], ptr[j + 1]):
                p *= T[i, cols[q]]
            s += p
        res[i] = s * scale
    return res


@njit(parallel=True, cache=True)
def apply_adj_csr(TT, ptr, cols, y, scale):
    n = ptr.size - 1
    m = y.size
    res = np.empty(n)
    for j in prange(n):
        s = 0.0
        a, b = ptr[j], ptr[j + 1]
        if a == b:                      # nu = 0: constant column
            for i in range(m):
                s += y[i]
        else:
            for i in range(m):
                p = y[i]
                for q in range(a, b):
                    p *= TT[cols[q], i]
                s += p
        res[j] = s * scale
    return res

def sparse_tensor_from_indices(J, Z, precomputed_table=chebyshev_table, normalization=None):
    """Build the univariate table and the sparse (CSR-like) description of the multi-indices."""
    Ja = np.asarray(J, dtype=np.int64)
    Za = np.asarray(Z, dtype=np.float64)
    m, d = Za.shape[0], Ja.shape[1]
    numax = Ja.max(axis=0)

    # Column offsets of each dimension in the (m, C) table. Degree 0 is never stored (T_0 = 1).
    offsets = np.concatenate(([0], np.cumsum(numax)))
    T = np.empty((m, offsets[-1]))
    for k in range(d):
        if numax[k] > 0:
            T[:, offsets[k]:offsets[k + 1]] = precomputed_table(Za[:, k], numax[k])[:, 1:]

    nz_rows, nz_dims = np.nonzero(Ja)                    # row-major order: grouped by multi-index
    cols = offsets[nz_dims] + Ja[nz_rows, nz_dims] - 1
    ptr = np.zeros(Ja.shape[0] + 1, dtype=np.int64)
    np.add.at(ptr, nz_rows + 1, 1)
    ptr = np.cumsum(ptr)

    if normalization is None:
        normalization = np.sqrt(m)
    return T, ptr, cols, 1.0 / normalization

@njit(parallel=True, cache=True)
def _columns(T, ptr, cols, U, scale):
    m = T.shape[0]
    res = np.empty((m, U.size))
    for i in prange(m):
        for jj in range(U.size):
            j = U[jj]
            p = scale
            for q in range(ptr[j], ptr[j + 1]):
                p *= T[i, cols[q]]
            res[i, jj] = p
    return res

def _npz(fname):
    return fname if fname.endswith(".npz") else fname + ".npz"
