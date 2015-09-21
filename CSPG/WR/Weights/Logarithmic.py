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

class Logarithmic:
    def __init__(self, c, alpha, stepsize=None):
        if stepsize is None:
            stepsize = 1

        self.name = 'Logarithmic {0}*j**{1}'.format(c, alpha)
        self.w    = c * np.repeat(np.arange(1,100000), stepsize)**alpha
