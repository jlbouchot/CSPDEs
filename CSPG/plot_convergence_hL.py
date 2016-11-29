from dolfin import *

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
import shelve


def plot_1D_hL(infile, outfile, L_to_plot, k_to_plot, which_to_plot = None):
    ### Load
    results = []
    for one_L in L_to_plot:
        cur_infile = infile + '_L{}'.format(one_L)
        print("Loading {0} ...".format(cur_infile))
        cur_results = sorted(shelve.open(cur_infile).values(), key=lambda r: r.L)
        results.append(cur_results[0]) # probably need a [0] here. 

    ### Plot
    sorted_results = sorted(results, key = lambda r: r.L)
    first_result = sorted_results[-1] # At that moment, first_result is a multi-level result
    #print len(first_result)
    d            = first_result.cspde_result[0].d # number of parameters
    spde_model   = first_result.spde_model
    epsilon      = first_result.epsilon


    d_plot       = len(k_to_plot) # k_to_plot corresponds to the number of pointwise convergences that we want to analyze
    results_plot = [r for r in results if r.L in L_to_plot]
    # Note: The number of grid points in the original grid can be recovered from spde_model.mesh_size - Obviously, every time we add a level, we multiply (every component) by two
    print len(results_plot)
    fig = plt.figure()
    # plt.suptitle("{0}: s={1}, epsilon={2}, {3} weights".format("Diffusion Model", s, epsilon, wr_model.weights.name))

    outer_grid = gridspec.GridSpec(int(np.ceil(d_plot/2.)), 2)

    l = []
    # First refine the mesh to a (much) finer level 
    spde_model.refine_mesh(2**(max(L_to_plot) + 2)) # Refine to 2 levels larger than the largest level we consider. !!!!!! To rechange back once the first tests work

    for (index, k) in enumerate(k_to_plot): 
        print("Plotting k={0} ...".format(k+1))

        # Reconstruct
        step = 0.02 # Change back to 0.02
        T = np.arange(-1+step, 1, step)
        Z = np.array([np.hstack((np.zeros(k), t, np.zeros(d-k-1))) for t in T])

		# Compute the solutions on that new mesh
        Y = [spde_model.samples(Z)] # This acts as the "True" solution
        for result in results_plot:
            # print result
			# Evaluation of the estimation / approximation in the ML case: 
            Y.append(result.wr_model.estimate_ML_samples(result.cspde_result, Z)) # Check that we actually have the level we are testing. 
			# Evaluation in the SL case: 
            # F_recon = result.wr_model.operator.create(np.array(result.cspde_result.J_s)[x != 0], Z, normalization=np.sqrt(result.cspde_result.m))
            # Y.append(F_recon.apply(x[x != 0]))
			# x       = result.cspde_result.result.x

        # Plot
        if which_to_plot is None:
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[index])

            ax = plt.Subplot(fig, inner_grid[0])
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
            for y in Y:
                l.append(ax.plot(T, y)[0])
            fig.add_subplot(ax)
            plt.title("$y_{0}$".format(k+1))

            ax = plt.Subplot(fig, inner_grid[1])
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
            for y in Y:
                ax.plot(T, Y[0]-y)
            fig.add_subplot(ax)
        else:
            inner_grid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_grid[index])
            ax = plt.Subplot(fig, inner_grid[0])
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

            if which_to_plot[index] == 'function':
                for y in Y:
                    l.append(ax.plot(T, y)[0])
                fig.add_subplot(ax)
                # plt.title("$z_{0}$".format(k+1))
            else:
                for y in Y:
                    ax.plot(T, Y[0]-y)
                fig.add_subplot(ax)
                # plt.title("$z_{0}$ difference".format(k+1))

    fig.legend(l, ["Truth"] + ["L={0}".format(int(round(result.L))) for result in results_plot],
               loc='lower center', bbox_to_anchor=(0.5, -0.05), bbox_transform=plt.gcf().transFigure,
               ncol=len(results_plot)+1, fancybox=True, shadow=True)


    all_axes = fig.get_axes()

    #show only the outside spines
    for ax in all_axes:
        if ax.is_last_row():
            plt.setp(ax.get_xticklabels(), visible=True)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)        

    if which_to_plot is None:
        plt.gcf().set_size_inches(17, int(4 * np.ceil(d_plot/2)))
    else:
        plt.gcf().set_size_inches(17, int(8 * np.ceil(d_plot/2)))

    # plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
