import numpy as np

## Add utilities to the path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'utilities'))
import utils as u

from .Result import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2015/09/21"
__lastmodified__ = "2026/07/21"

def primaldual(Operator, P_Fstar, P_G, theta, tau, sigma, E, eta, maxiter):
    # Initialize variables
    xi    = np.zeros(Operator.m)
    x_new = np.zeros(Operator.n)
    x_old = np.zeros(Operator.n)

    k = 0

    t = u.time_things()

    # Compute solution
    while E(x_new, xi) > eta:
        x_new = P_G(tau, x_old - tau*Operator.apply_adj(xi)).T
        xi    = P_Fstar(sigma, xi + sigma*Operator.apply(x_new + theta*(x_new - x_old)))
        x_old = x_new

        k += 1

        if k > maxiter:
            print('Primal-Dual did not converge after {0} steps.'.format(k))
            break

    t = u.time_things(t)

    return Result(x_new, k, 'Primal Dual', k <= maxiter, t[0], t[1], t[2])
