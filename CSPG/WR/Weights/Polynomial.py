import numpy as np

class Polynomial:
    def __init__(self, c, alpha, stepsize=None):
        if stepsize is None:
            stepsize = 1

        if c == 1:
            self.name = 'Polynomial $j^{{{0}}}$'.format(alpha)
        else:
            self.name = 'Polynomial ${0} \cdot j^{{{1}}}$'.format(c, alpha)
            
        self.w    = c * np.repeat(np.arange(1,100000), stepsize)**alpha
