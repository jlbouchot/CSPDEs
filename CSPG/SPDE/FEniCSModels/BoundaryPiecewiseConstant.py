from dolfin import *

import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2017/11/15"

class BoundaryPiecewiseConstant:
    def __init__(self, d, local_var = 0.5, c = None): # c corresponds to the average (in expectation) Dirichlet BC.
        self.num_params = 2*d # this is used only for 2D problems at the moment with d pieces on the left boundary and d pieces on the right one.

        if c is None:
            self.c = 5
        else:
            self.c = c

        self.local_var = local_var

    def __call__(self, x, z): # we could have some more complicated things by having some x-variations and similar
        print "Calling Boundary Pieccewise Coefficient !!!!!"
        return [self.c]*self.num_params + local_var*z # Return a list of the local parameters
