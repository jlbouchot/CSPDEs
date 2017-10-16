from dolfin import *

import numpy as np
from scipy.special import zeta

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2017/10/16"

class CosineCoef2D:
    def __init__(self, d, alpha, c = None):
        self.num_params = d

        if c is None:
            self.c = zeta(alpha, 1) * 2 + 0.1
        else:
            self.c = c

        self.alpha = alpha

    def __call__(self, x, z):
        r = 0
        for k in range(int(len(z))):
            r +=     z[k]   * cos(np.pi * (k+1) * (x[0]**2 + x[1]**2 + x[2]**2)) / (k+1)**self.alpha
        return self.c + r
