from dolfin import *

import numpy as np
from scipy.special import zeta


class TrigCoefficient:
    def __init__(self, d, alpha, c = None):
        self.num_params = 2*d

        if c is None:
            self.c = zeta(alpha, 1) * 2 + 0.1
        else:
            self.c = c

        self.alpha = alpha

    def __call__(self, x, z):
        r = 0

        for k in range(int(len(z)/2)):
            r +=     z[2*k]   * cos(np.pi * (k+1) * x[0]) / (k+1)**self.alpha

            if 2*k+1 < len(z):
                r += z[2*k+1] * sin(np.pi * (k+1) * x[0]) / (k+1)**self.alpha

        return self.c + r
