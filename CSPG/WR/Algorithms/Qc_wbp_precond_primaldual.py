import numpy as np

from Result     import *
from primaldual import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

class Qc_wbp_precond_primaldual:
    def __init__(self, theta, eta_pd, maxiter):
        self.theta   = theta
        self.eta_pd  = eta_pd
        self.maxiter = maxiter
    

    def __call__(self, Operator, y, w, s, eta):
        ### Error function for indication whether an approximation is good enough yet
        def E(x, xi):
            return np.linalg.norm(Operator.apply(x) - y)

        ### Functions related to F
        def P_Fstar(sigma, x):
            if np.linalg.norm(x/sigma - y) <= eta:
                return 0

            return (1 - eta/np.linalg.norm(x/sigma - y)) * (x - sigma*y);

        ### Functions related to G
        def P_G(tau, x):
            r = np.zeros(len(x))

            for k in range(len(x)):
                r[k] = P_abs(tau[k]*w[k], x[k])

            return r

        def P_abs(tau, z):
            if abs(z) >= tau:
                return np.sign(z) * (abs(z) - tau)

            return 0


        sigma = 1./np.array([np.sum(np.abs(Operator.A[:, k])) for k in range(Operator.n)])
        tau   = 1./np.array([np.sum(np.abs(Operator.A[l, :])) for l in range(Operator.m)])

        # sigma = 1./np.array([np.sum(np.abs(Operator.apply(np.eye(Operator.n, 1, -k))))     for k in range(Operator.n)])
        # tau   = 1./np.array([np.sum(np.abs(Operator.apply_adj(np.eye(Operator.m, 1, -l)))) for l in range(Operator.m)])

        print "Done computing preconditioning matrices, now actually solving ..."
        r = primaldual(Operator, P_Fstar, P_G, self.theta, sigma, tau, E, self.eta_pd, self.maxiter)

        return Result(r.x, r.iterations, 'Preconditioned Quadratically Constrained ' + r.methodname)
