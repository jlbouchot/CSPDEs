import numpy as np
from numba import njit,prange
__author__ = ["Falk Pulsmeyer", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/02/22"

class operator_from_matrix_Alt:
    def __init__(self, model, A, univariate_A,J):
        self.model = model
        self.A     = A
        self.univariate_A = univariate_A
        self.m     = self.univariate_A[0].shape[0]
        self.n     = len(J)
        self.J     = np.array(J)

    def apply(self, x):
        t=self.apply_alt(x)
        return t


    def apply_alt(self, x):
        res=np.empty(self.univariate_A[0].shape[0])
        t = apply_alt_numba(self.univariate_A,res,x,self.J)
        return t


    def apply_adj(self, x):
        t=self.apply_adj_alt(x)
        return t

    def apply_adj_alt(self, y):
        res = np.empty(len(self.J)) #number of the samples
        t = apply_adj_alt_numba(self.univariate_A,res,y,self.J)
        return t

    def genSubMatrix(self, U):
        t=self.genSubMatrix_Alt(U)
        return t

    def genSubMatrix_Alt(self, U):
        res=np.zeros((self.univariate_A[0].shape[0], len(U)))
        t = genSubMatrix_Alt_Numba(self.univariate_A,res,U,self.J)
        return t

    def save(self, data_mtx):
        np.save(data_mtx, [self.A,self.univariate_A,self.J])



def matrix_from_tensor_indices(J, Z, base, normalization=None):
    Za, Ja = np.array(Z), np.array(J)

    # A = np.array([np.array([np.prod(base(z_row, j)) for j in Ja]) for z_row in Za])
    # A = np.reshape(np.fromiter((np.prod(base(z_row, j)) for j in Ja for z_row in Za), np.float), (len(Za), len(Ja)), len(Za)*len(Ja))
    A = np.zeros((3,3))
    if normalization is None:
        normalization = np.sqrt(np.size(A, 0))

    return A/normalization


def univ_tensor_from_tensor_indices(J, Z, base, normalization=None):
    Za, Ja = np.array(Z), np.array(J)
    numax=np.amax(Ja,0) #amax returns the maximum entry accodirng to the specified axis in this case the 0-th one

    univariate_A = []
    if normalization is None:
        normalization = np.sqrt(len(Za))
    normalization = normalization**(1/len(numax))
    for i in range(0,len(numax)):#dimension
        tmp=np.empty((len(Za), numax[i]+1))
        for k in range(0,numax[i]+1): #degree
            for j in range(0,len(Za)):#sample
                tmp[j,k] = base(Za[j,i],k)/normalization
        univariate_A.append(tmp)
    return univariate_A

@njit
def apply_alt_numba(wert, res, x,J):
    for i in prange(wert[0].shape[0]):
        my_sum = 0
        for j in prange(len(x)):
            my_prod=1
            for k in range(0, len(wert)):
                my_prod*=wert[k][i][J[j][k]]
            my_sum+=my_prod*x[j]
        res[i]=my_sum
    # return res/(np.sqrt(wert[0].shape[0]))
    return res

@njit
def apply_adj_alt_numba(wert, res, y, J):
    for i in range(0, len(J)):
        my_sum = 0
        for j in range(0, wert[0].shape[0]):
            my_prod = 1
            for k in range(0, len(wert)):
                my_prod *= wert[k][j][J[i][k]]
            my_sum += my_prod * y[j]
        res[i] = my_sum
    # return res/(np.sqrt(wert[0].shape[0]))
    return res

@njit
def genSubMatrix_Alt_Numba(wert,res,U,J):
    for i in range(0,res.shape[0]):
        for j in range(0,len(U)):
            tmp = 1
            for k in range(0, len(wert)):
                tmp*=wert[k][i][J[U[j]][k]]
            res[i][j] = tmp
    # return res/(np.sqrt(wert[0].shape[0]))
    return res
