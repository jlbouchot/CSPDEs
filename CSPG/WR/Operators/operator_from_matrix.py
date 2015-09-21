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

class operator_from_matrix:
    def __init__(self, model, A):
        self.model = model
        self.A     = A
        self.m     = np.size(A, 0)
        self.n     = np.size(A, 1)

    def apply(self, x):
        return np.dot(self.A, x)

    def apply_adj(self, x):
        return np.dot(self.A.T, x)

    def genSubMatrix(self, U):
        return self.A[:, U]


def matrix_from_tensor_indices(J, Z, base, normalization=None):
    Za, Ja = np.array(Z), np.array(J)

    # A = np.array([np.array([np.prod(base(z_row, j)) for j in Ja]) for z_row in Za])
    A = np.reshape(np.fromiter((np.prod(base(z_row, j)) for j in Ja for z_row in Za), np.float), (len(Za), len(Ja)), len(Za)*len(Ja))

    if normalization is None:
        normalization = np.sqrt(np.size(A, 0))

    return A/normalization
