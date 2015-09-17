import numpy as np

class Logarithmic:
    def __init__(self, c, alpha, stepsize=None):
        if stepsize is None:
            stepsize = 1

        self.name = 'Logarithmic {0}*j**{1}'.format(c, alpha)
        self.w    = c * np.repeat(np.arange(1,100000), stepsize)**alpha
