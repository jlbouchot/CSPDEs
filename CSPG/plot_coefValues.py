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

from collections import namedtuple

MCResults = namedtuple('MCResults', ['nu', 'MC_est', 'l1_est'])


def plot_coef_values(infile, outfile, ratio = 1./10.): #, L_to_plot, v_to_plot, k_to_plot, which_to_plot = None):
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






def plot_coef_from_file(infile, imgfile):
    data = shelve.open(infile).values() # This contains an array of an MCResults structure with fields .nu (an np.array), MC_est (containing the MC estimations!), l1_est (the estimation resulting from CS)

    
    eps_not_zero = 1e-22 # make sure we can take a log!
    # true_coefs = log10(true_coefs + eps_not_zero)
    true_coefs = [math.log10(np.abs(dummy.MC_est) + eps_not_zero) for dummy in data]
#    computed_coefs = log10(computed_coefs + eps_not_zero)
    computed_coefs = [math.log10(np.abs(dummy.l1_est) + eps_not_zero) for dummy in data]

    nb_coefs = len(true_coefs)
    
    # Make sure the two vectors to plot have the same ordering
    sorting_permutations = sorted(range(len(true_coefs)), key=lambda k: true_coefs[k])
    sorted_computed_coefs = [computed_coefs[i] for i in sorting_permutations]
    sorted_true_coefs = [true_coefs[i] for i in sorting_permutations]

    nb_to_display = nb_coefs
    # red dashes, blue squares and green triangles
    plt.plot(range(0,nb_to_display), sorted_true_coefs, color='r', marker='^', label='MC Estimations')
    plt.plot(range(0,nb_to_display), sorted_computed_coefs, color='b', marker='s', label='Multi-level Estimations')
    # handler_true_coefs = plt.plot(range(0,nb_to_display), sorted_true_coefs[sorting_idx[0:nb_to_display]], color='r', marker='^', label='True coefficients')
    # handler_estim_coefs = plt.plot(range(0,nb_to_display), sorted_computed_coefs[sorting_idx[0:nb_to_display]], color='b', marker='s', label='Multi-level estimations')
    plt.xlabel('Coefficient index')
    plt.ylabel('Coefficients log-magnitudes')
    # plt.legend(handles = [handler_true_coefs, handler_estim_coefs])
    plt.legend()
    
    plt.savefig(imgfile, bbox_inches="tight")
    plt.show()





def plot_coef_from_file2(infile, imgfile):
    data = sorted(shelve.open(infile).values(), key = lambda r: np.abs(r.MC_est), reverse = True ) # This contains an array of an MCResults structure with fields .nu (an np.array), MC_est (containing the MC estimations!), l1_est (the estimation resulting from CS)

    
    eps_not_zero = 1e-22 # make sure we can take a log!
    # true_coefs = log10(true_coefs + eps_not_zero)
    true_coefs = [math.log10(np.abs(dummy.MC_est) + eps_not_zero) for dummy in data]
#    computed_coefs = log10(computed_coefs + eps_not_zero)
    computed_coefs = [math.log10(np.abs(dummy.l1_est) + eps_not_zero) for dummy in data]

    nb_coefs = len(true_coefs)
    
    # Make sure the two vectors to plot have the same ordering
    sorting_permutations = sorted(range(len(true_coefs)), key=lambda k: true_coefs[k])
    sorted_computed_coefs = [true_coefs[i] for i in sorting_permutations]
    sorted_true_coefs = [computed_coefs[i] for i in sorting_permutations]

    nb_to_display = nb_coefs
    # red dashes, blue squares and green triangles
    plt.plot(range(0,nb_to_display), sorted_true_coefs, color='r', marker='^', label='MC Estimations')
    plt.plot(range(0,nb_to_display), sorted_computed_coefs, color='b', marker='s', label='Multi-level Estimations')
    # handler_true_coefs = plt.plot(range(0,nb_to_display), sorted_true_coefs[sorting_idx[0:nb_to_display]], color='r', marker='^', label='True coefficients')
    # handler_estim_coefs = plt.plot(range(0,nb_to_display), sorted_computed_coefs[sorting_idx[0:nb_to_display]], color='b', marker='s', label='Multi-level estimations')
    plt.xlabel('Coefficient index')
    plt.ylabel('Coefficients log-magnitudes')
    # plt.legend(handles = [handler_true_coefs, handler_estim_coefs])
    plt.legend()
    
    plt.savefig(imgfile, bbox_inches="tight")
    plt.show()





# Using the MCestimations from smallDim6_PWLD_whtp_var_big_grid2000_L4_g105_datconstant20_MCResults_1500000
# plot with smallDim6_PWLD_womp_var_small_grid2000_L3_g105_datconstant20
# or smallDim6_PWLD_womp_var_small_grid2000_L3_g1055_datconstant20
def better_plot_coefs(MC_file, l1_file, imgfile): 
    data = sorted(shelve.open(MC_file).values(), key = lambda dummy: np.abs(dummy.MC_est), reverse = True)
    
    eps_not_zero = 1e-22 # make sure we can take a log!
    MC_values = [math.log10(np.abs(dummy.MC_est) + eps_not_zero) for dummy in data]

    # # The dictionary based representation i redundant, but it allows for a somewhat better readability as well as more adaptability in other formats of results
    # dict_MC = dict(zip([dummy.nu for dummy in data], [dummy.MC_est for dummy in data]))


    l1_estimations, J = gcc(l1_file, [], J = None)
    l2_estimations, J_l2 = gtc(l1_file, None, J)
    # dict_l1 = {}
    l1_values = []
    l2_values = []
    for a_nu in [dummy.nu for dummy in data]:
        find_idx = np.array(np.sum(J == a_nu, axis = 1))
        cur_idx_estimations = find_idx.tolist().index(len(a_nu))
        # dict_l1[a_nu] = l1_esimations[cur_idx_estimations]
        l1_values.append(math.log10(np.abs(l1_estimations[cur_idx_estimations]) + eps_not_zero))

        find_idx = np.array(np.sum(J_l2 == a_nu, axis = 1))
        cur_idx_estimations = find_idx.tolist().index(len(a_nu))
        # dict_l1[a_nu] = l1_esimations[cur_idx_estimations]
        l2_values.append(math.log10(np.abs(l2_estimations[cur_idx_estimations]) + eps_not_zero))

    nb_coefs = len(MC_values)
    nb_to_display = nb_coefs
    # red dashes, blue squares and green triangles
    plt.plot(range(0,nb_to_display), MC_values, color='r', marker='^', label='MC Estimations')
    plt.plot(range(0,nb_to_display), l1_values, color='b', marker='s', label='Multi-level Estimations')
    plt.plot(range(0,nb_to_display), l2_values, color='y', marker='x', label='L2 Estimations')
    # plt.plot(range(0,nb_to_display), dict_MC.values(), color='r', marker='^', label='MC Estimations')
    # plt.plot(range(0,nb_to_display), dict_l1.values(), color='b', marker='s', label='Multi-level Estimations')
    # handler_true_coefs = plt.plot(range(0,nb_to_display), sorted_true_coefs[sorting_idx[0:nb_to_display]], color='r', marker='^', label='True coefficients')
    # handler_estim_coefs = plt.plot(range(0,nb_to_display), sorted_computed_coefs[sorting_idx[0:nb_to_display]], color='b', marker='s', label='Multi-level estimations')
    plt.xlabel('Coefficient index')
    plt.ylabel('Coefficients log-magnitudes')
    # plt.legend(handles = [handler_true_coefs, handler_estim_coefs])
    plt.legend()
    
    plt.savefig(imgfile, bbox_inches="tight")
    plt.show()
