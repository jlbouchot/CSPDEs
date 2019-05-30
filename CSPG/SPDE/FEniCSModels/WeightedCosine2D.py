from dolfin import *

import numpy as np
from scipy.special import zeta

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institut of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/05/30"

class WeightedCosine2D:


    def __init__(self, d, alpha, imp = 1, c = None, weight_cst = 0.5):
        """
        Parameters
        ----------
        d : int
            The number of Cosines in the expansion
        alpha : float
            the power of the decay of the cosine
        imp : float, optional
            The (global) importance of the fluctuations with respect to the main field (Default = 1)
        c : float, optional
            The value of the (Constant) mean field. (Default depends on the alpha via zeta function)
        weight_cst : float (positive), optional
            The amount to which the current fluctuation is taken into account (Default = 1/2)
        """

        self.num_params = d

        if c is None:
            self.c = zeta(alpha, 1) * 2 + 0.1
        else:
            self.c = c

        self.alpha = alpha
        self.imp = imp
        self.weight_cst = weight_cst

    def __call__(self, x, z):
        r = 0
        for k in range(int(len(z))):
            r +=     z[k]   * cos(np.pi * (k+1) * (x[0]**2+x[1]**2 )**0.5) / (self.weight_cst*(k+1)+1)**self.alpha
        return self.c + self.imp*r
