import numpy as np

from Result import *


def primaldual(Operator, P_Fstar, P_G, theta, tau, sigma, E, eta, maxiter):
    # Initialize variables
    xi    = np.zeros(Operator.m)
    x_new = np.zeros(Operator.n)
    x_old = np.zeros(Operator.n)

    k = 0

    # Compute solution
    while E(x_new, xi) > eta:
        x_new = P_G(tau, x_old - tau*Operator.apply_adj(xi)).T
        xi    = P_Fstar(sigma, xi + sigma*Operator.apply(x_new + theta*(x_new - x_old)))
        x_old = x_new

        k += 1

        if k >= maxiter:
            print('Primal-Dual did not converge after {0} steps.'.format(k))
            break

    return Result(x_new, k, 'Primal Dual')
