import numpy as np

def weighted_quasi_abslargest(x, s, w):
    # Get indices of elements sorted in descending order
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
