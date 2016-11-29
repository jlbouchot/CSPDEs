import numpy as np

from dolfin import *

import sys
import shelve
import datetime
import time

import os.path

import scipy.integrate as spi

from compare_true_coefs import get_computed_coefs as gcc
from iterative_solution import compute_true_avg_alternate as cta
from multivariate_chebpoly import multi_var_chebpoly as mvc


# sp.integrate.quad(f,a,b)
# cta(ys, rhs, variability, mean_field)

# first get the computed coefficients as well as the set of all potential candidate coefficients
computed_coefs, other_J = gcc(infile, None, J) # other_J should be the same as J.

# For every multi index in J,
# Compute (numerical integration) the estimated coefficient
# Compare! 
c_nus = np.zeros(len(other_J))
for idx, one_nu in enumerate(other_J):
    # Construct the current Chebychef polynomial with the associated weight
    for nu_j in one_nu:
        if nu_j

spi.nquad(lambda x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12: cta(np.array([x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]), 1.0, 1.0/13, 5.0 )*mvc([x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12],a_nu),[[-1,1]]*13)
# compute_true_avg(ys, rhs, variability, mean_field):

spi.nquad(lambda x0, x1: cta(np.array([x0,x1,0.01,-0.1,-0.02,0,0,0.02,0,0,0,0,0]).transpose(), 1.0, 1.0/13, 5.0 )*mvc([x0,x1,0.01,-0.1,-0.02,0,0,0.02,0,0,0,0,0],a_nu),[[-1,1]]*2)


import numpy
import scipy.integrate
import math

def w(r, theta, phi, alpha, beta, gamma):
    return(-math.log(theta * beta))

def integrand(phi, alpha, gamma, r, theta, beta):
    ww = w(r, theta, phi, alpha, beta, gamma)
    k = 1.
    T = 1.
    return (math.exp(-ww/(k*T)) - 1.)*r*r*math.sin(beta)*math.sin(theta)

# limits of integration

def zero(x, y=0):
    return 0.

def one(x, y=0):
    return 1.

def pi(x, y=0):
    return math.pi

def twopi(x, y=0):
    return 2.*math.pi

# integrate over phi [0, Pi), alpha [0, 2 Pi), gamma [0, 2 Pi)
def secondIntegrals(r, theta, beta):
    res, err = scipy.integrate.tplquad(integrand, 0., 2.*math.pi, zero, twopi, zero, pi, args=(r, theta, beta))
    return res

# integrate over r [0, 1), beta [0, 2 Pi), theta [0, 2 Pi)
def integral():
    return scipy.integrate.tplquad(secondIntegrals, 0., 2.*math.pi, zero, twopi, zero, one)
