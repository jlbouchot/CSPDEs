from dolfin import *
import ufl

import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/22"

class PartitionUnityConstantCoefficient:
    def __init__(self, d = None, abar = None, coefs = None, parts = None):
        if abar is None: 
            self.abar = 1
        else:
            self.abar = abar
        if d is None and coefs is None and parts is None:
            self.num_params = 10
            self.coefs = np.absolute(abar)*1.0/(self.num_params+1)*np.ones(self.num_params)
            self.parts = np.linspace(1,d,d)*1.0/self.num_params
        # At least one of the parameters is defined
        if d is not None: 
            self.num_params = d
            if parts is None:
                self.parts = np.linspace(1,d,d)*1.0/self.num_params
                if coefs is None:
                    self.coefs = np.absolute(abar)*1.0/(self.num_params+1)*np.ones(self.num_params)
                else: 
                    self.coefs = coefs
            else:
                self.parts = parts 
                if coefs is None:
                    self.coefs = np.absolute(abar)*1.0/	(self.num_params+1)*np.ones(self.num_params)
                else: 
                    self.coefs = coefs
        # We should now check that the user hasn't made a mistake when passing argument
        # We should also check that the last element of the parts is 1. 
        
    def __call__(self, x, z):

        def local_variation(k, x):
            print x[0]
            print self.parts[k]
            return conditional(le(ufl.geometry.evaluate(x), self.parts[k]), self.coefs[k], local_variation(k+1,x))
            # if x <= self.parts[k]:
            #     return self.coefs[k]
            # else:
            #     return local_variation(k+1,x)

        # print self.coefs
        # print self.abar
        # print self.parts
        # print x
        # print type()
        return self.abar + local_variation(0, x)