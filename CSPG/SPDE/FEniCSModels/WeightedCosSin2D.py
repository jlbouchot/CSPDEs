from dolfin import *

import numpy as np

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015 - 2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/30"
__lastmodified__ = "2026/10/30"

class WeightedCosSin2D:
    '''
    a(x, z) = c + imp * sum_{k=1}^{d/2} [z_{2k-1} cos(pi k r) + z_{2k} sin(pi k r)] / (weight_cst k + 1)^alpha,  r = |x|

    Same constructor as WeightedCosine2D, but d is the number of *parameters* (two per frequency, d even).
    Both functions of a frequency have the same decay, so the duplicated weights of the drivers
    (v = gamma * [k for k in 1..d/2 for _ in (0,1)]^exponent) are the right ones for this coefficient.
    '''

    def __init__(self, d, alpha, imp = 1, c = None, weight_cst = 0.5):
        assert d % 2 == 0, "WeightedCosSin2D needs an even number of parameters"
        self.num_params = d
        self.c = 2.0 if c is None else c
        self.alpha = alpha
        self.imp = imp
        self.weight_cst = weight_cst

    def __call__(self, x, z):
        r = 0
        rad = (x[0]**2 + x[1]**2)**0.5
        for j in range(int(len(z))):
            k = j // 2 + 1
            f = cos if j % 2 == 0 else sin
            r += z[j] * f(np.pi * k * rad) / (self.weight_cst*k + 1)**self.alpha
        return self.c + self.imp*r
