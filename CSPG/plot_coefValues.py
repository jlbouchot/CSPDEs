from dolfin import *

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from compare_true_coefs import get_true_coefs as gtc
from compare_true_coefs import get_computed_coefs as gcc

import sys
import shelve

import math


def plot_coef_values(infile, outfile, ratio = 1): #, L_to_plot, v_to_plot, k_to_plot, which_to_plot = None):
    ### Load
    true_coefs, J = gtc(infile, None, None)
    computed_coefs, other_J = gcc(infile, None, J) # other_J should be the same as J.
    
    J = np.asarray(J)
    other_J = np.asarray(other_J)
    d = len(J[0,:]) # number of parameters -- make sure that's the way to go!

    eps_not_zero = 1e-22 # make sure we can take a log!
    # true_coefs = log10(true_coefs + eps_not_zero)
    true_coefs = [math.log10(np.abs(a_coefs) + eps_not_zero) for a_coefs in true_coefs]
#    computed_coefs = log10(computed_coefs + eps_not_zero)
    computed_coefs = [math.log10(np.abs(a_coef) + eps_not_zero) for a_coef in computed_coefs]

    nb_coefs = len(true_coefs)
    
    # Make sure the two vectors to plot have the same ordering
    sorted_computed_coefs = np.zeros(nb_coefs)
    sorted_true_coefs = np.zeros(nb_coefs)
    for idx, cur_nu in enumerate(J): #xrange(0,d):
        find_idx = np.array(np.sum(J == cur_nu, axis = 1))
        cur_idx_true = find_idx.tolist().index(d)
        find_idx = np.array(np.sum(other_J == cur_nu, axis = 1))
        cur_idx_estimations = find_idx.tolist().index(d)
        sorted_computed_coefs[cur_idx_estimations] = computed_coefs[cur_idx_estimations]
        sorted_true_coefs[cur_idx_true] = true_coefs[cur_idx_true]
    
    # # Now sort and keep the permutation somewhere
    sorting_idx = np.argsort(np.abs(sorted_true_coefs))
    # plt.plot([1,2,3,4], [1,4,9,16], 'ro')
    # plt.axis([0, 6, 0, 20])
    # plt.show()

    # or
    nb_to_display = int(np.floor(ratio*nb_coefs))
    # red dashes, blue squares and green triangles
    plt.plot(range(0,nb_to_display), sorted_true_coefs[sorting_idx[0:nb_to_display]], color='r', marker='^', label='True coefficients')
    plt.plot(range(0,nb_to_display), sorted_computed_coefs[sorting_idx[0:nb_to_display]], color='b', marker='s', label='Multi-level estimations')
    # handler_true_coefs = plt.plot(range(0,nb_to_display), sorted_true_coefs[sorting_idx[0:nb_to_display]], color='r', marker='^', label='True coefficients')
    # handler_estim_coefs = plt.plot(range(0,nb_to_display), sorted_computed_coefs[sorting_idx[0:nb_to_display]], color='b', marker='s', label='Multi-level estimations')
    plt.xlabel('Coefficient index')
    plt.ylabel('Coefficients log-magnitudes')
    # plt.legend(handles = [handler_true_coefs, handler_estim_coefs])
    plt.legend()
    
    plt.savefig(outfile, bbox_inches="tight")
    plt.show()
