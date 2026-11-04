import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/10/28"

def weighted_quasi_abslargest(x, s, w):
    '''
    Find the approximate indices of the elements with largest weighted absolute values
        x: input vector
        s: target weighted sparsity
        w: weight sequence
    '''
    sortIndex = np.argsort((np.abs(x)) * (w**(-1)))[::-1]

    k      = 0
    w_test = 0
    w_cur  = 0

    while True:
        w_test = w_cur + w[sortIndex[k]]**2

        if s < w_test:
            break

        w_cur = w_test
        k     = k + 1


    i    = sortIndex[0:k]

    r    = np.zeros_like(x)
    r[i] = x[i]

    return r, i

def weighted_quasi_abslargest_fast(x, s, w):
    '''
    Vectorised weighted thresholding, drop-in replacement for weighted_quasi_abslargest

    only the candidates that can possibly be selected are sorted: since the weights
    satisfy w_j >= w_min, at most floor(s / w_min^2) entries fit in the budget, so an
    O(N) argpartition followed by a sort of those few entries replaces the O(N log N)
    full argsort on vectors of length N ~ 1e5-1e6;
    '''
    score = np.abs(x) / w
    n = score.size
    kmax = int(min(n, np.floor(s / np.min(w) ** 2) + 1))

    if kmax < n:
        cand = np.argpartition(-score, kmax - 1)[:kmax]
        order = cand[np.argsort(score[cand])[::-1]]
    else:
        order = np.argsort(score)[::-1]

    k = int(np.searchsorted(np.cumsum(w[order] ** 2), s, side="right"))
    i = order[:k]

    r = np.zeros_like(x)
    r[i] = x[i]
    return r, i