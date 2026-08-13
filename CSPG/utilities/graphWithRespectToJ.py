# Graph total compute time, PDE only compuyte time, l1recovery only compute time with respect to J

from dolfin import *


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pandas as pd

import sys
import shelve

import os.path 
# thanks python >= 3.5 for the next thing! 
import pathlib # This allows to avoid the if path.exists -> mkdir thing. Hurray!

from collections import namedtuple

## Add utilities to the path
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import Check_ML
from utilities.utils import getComputeTimes, getGroundTruthFromModel

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/08/13"

# Keep in mind the results' folders will all have the following form: 
# TODO: Add proper naming conventions

# Remember the things? Surely there is another way to do this!
TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery', 't_J'])

nbDim = 2
target_mesh_size = [3000]*nbDim
if "DEBUG_MODE" in os.environ:
	Js_to_display = [1,2,3] # [1,2,3,4,5] # Note that J = 5 corresponds to the SL appraoch
else:
	Js_to_display = [1,2,3,4,5] # Note that J = 5 corresponds to the SL appraoch

if "DEBUG_MODE" in os.environ:
	core_folder_name = os.path.join('UMFPACKresults', 'Exp3_debug')
else:
	core_folder_name = os.path.join('results', 'Exp3_debug')
base_fname = 'J_val_'
# core_folder_name = 'Exp4H020Dim2WCosine10InfluenceTarget5J'
# # core_folder_name = 'Exp4H020Dim2WCosine10InfluenceJ'
# fname_to_read = 'WeightedCosine2D' # This is an unhappy mistake in my code which makes all file to have the same name. Luckily, They are all saved in separate folders. 
cfg_fname = 'config_file.txt' # This contains all the details from the experiments. I don't think we need it for graphing, but who knows. 

## Placeholder macros
## For some reasons, I must have missed a python update or something, but I can't access the various locations of the fields of the classes I saved. 
##
## TR for the TestResult macros
#TR_SPDE_MODEL_POS = 0
#TR_WR_MODEL_POS = 1
#TR_EPSILON_POS = 2
#TR_L_POS = 3
#TR_CR_POS = 4
## CR for CSPDEResult macros
#CR_J_S_POS = 0
#CR_N_POS = 1
#CR_S_POS = 2
#CR_M_POS = 3
#CR_D_POS = 4
#CR_Z_POS = 5
#CR_Y_POS = 6
#CR_A_POS = 7
#CR_W_POS = 8
#CR_RESULT_POS = 9
#CR_T_SAMPLES_POS = 10
#CR_T_MATRIX_POS = 11
#CR_T_RECOVERY_POS = 12


#############
## Check if some ground truth data exist
#############
path_to_GT = 'groundTruth'

# Load and plot results for all d's, one after the other 
results_all = {}
for idx, oneJ in enumerate(Js_to_display): 
	cur_path_to_file = core_folder_name
	fname_to_read = base_fname + str(Js_to_display[idx])
	print("Loading {0} from folder {1} ...".format(fname_to_read, cur_path_to_file))
	cur_results = sorted(shelve.open(os.path.join(cur_path_to_file,fname_to_read)).values(), key=lambda r: r.L)
	if len(cur_results) > 0: 
		# results_all.append(cur_results[0]) # This is the Check_ML.TestResult tuple
		results_all[oneJ] = cur_results[0]
	else: 
		# results_all.append(cur_results)
		results_all[oneJ] = cur_results

finest_result	= results_all[Js_to_display[0]]
d 		= finest_result.cspde_result[0].d # number of parameters
nb_tests 	= 1000 # This should be sufficient
wr_model 	= finest_result.wr_model
spde_model 	= finest_result.spde_model

spde_model.set_mesh_size(target_mesh_size)
print("Getting ground truth from model with width {}".format(spde_model.mesh_size[0]))
y_GT, Z = getGroundTruthFromModel(spde_model, wr_model, d, nb_tests, GTfolder = 'GTresults_d' + str(d) + 'H' + str(target_mesh_size[0]), GTfilenames = 'GT_')


l2error = {} 
linferror = {} 
computeTimePDEWC = {} 
computeTimeRecoveryWC = {} 
computeTimeMtxWC = {} 
computeTimeJsWC = {} 
computeTotalTimeWC = {} 

computeTimePDEUser = {} 
computeTimeRecoveryUser = {} 
computeTimeMtxUser = {} 
computeTimeJsUser = {} 
computeTotalTimeUser = {} 

computeTimePDESys = {} 
computeTimeRecoverySys = {} 
computeTimeMtxSys = {} 
computeTimeJsSys = {} 
computeTotalTimeSys = {}

for oneJ in Js_to_display: 
	print("Computing results for J = {}".format(oneJ))
	# spde_model = results_all[oneJ].spde_model 
	# spde_model.set_mesh_size(target_mesh_size)
	# y_estimated = wr_model.estimate_ML_samples(first_result[0].cspde_result, Z)
	y_estimated = results_all[oneJ].wr_model.estimate_ML_samples(results_all[oneJ].cspde_result, Z)
	l2error[oneJ] = np.linalg.norm(y_estimated - y_GT)
	linferror[oneJ] = np.linalg.norm(y_estimated - y_GT, ord=np.inf)
	curPdeTimes, curRecoveryTimes, curMtxTimes, curJsTimes = getComputeTimes(results_all[oneJ].cspde_result)
	computeTotalTimeWC[oneJ] = curPdeTimes[0] + curRecoveryTimes[0] + curMtxTimes[0] + curJsTimes[0]
	computeTimePDEWC[oneJ] = curPdeTimes[0] / computeTotalTimeWC[oneJ]
	computeTimeRecoveryWC[oneJ] = curRecoveryTimes[0] / computeTotalTimeWC[oneJ]
	computeTimeMtxWC[oneJ] = curMtxTimes[0] / computeTotalTimeWC[oneJ]
	computeTimeJsWC[oneJ] = curJsTimes[0] / computeTotalTimeWC[oneJ]

	computeTotalTimeUser[oneJ] = curPdeTimes[1] + curRecoveryTimes[1] + curMtxTimes[1] + curJsTimes[1]
	computeTimePDEUser[oneJ] = curPdeTimes[1] / computeTotalTimeUser[oneJ]
	computeTimeRecoveryUser[oneJ] = curRecoveryTimes[1] / computeTotalTimeUser[oneJ]
	computeTimeMtxUser[oneJ] = curMtxTimes[1] / computeTotalTimeUser[oneJ]
	computeTimeJsUser[oneJ] = curJsTimes[1] / computeTotalTimeUser[oneJ]

	computeTotalTimeSys[oneJ] = curPdeTimes[2] + curRecoveryTimes[2] + curMtxTimes[2] + curJsTimes[2]
	computeTimePDESys[oneJ] = curPdeTimes[2] / computeTotalTimeSys[oneJ]
	computeTimeRecoverySys[oneJ] = curRecoveryTimes[2] / computeTotalTimeSys[oneJ]
	computeTimeMtxSys[oneJ] = curMtxTimes[2] / computeTotalTimeSys[oneJ]
	computeTimeJsSys[oneJ] = curJsTimes[2] / computeTotalTimeSys[oneJ]


colours = plt.cm.rainbow(np.linspace(0, 1, len(Js_to_display)+1))
markers = ["o", "v", "^", "<", ">", "s", "8"]


# plt.figure()
# cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
# scatter = []
# for idx,oneJ in enumerate(Js_to_display): 
# 	scatter.append(plt.scatter(np.log10(computeTotalTime[oneJ]), np.log10(linferror[oneJ]), color = colours[idx], marker=markers[idx] ))
# plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
# plt.xlabel('Total compute time ($log_{10}$ scale)')
# classes = ["J = " + str(oneJ) for oneJ in Js_to_display]
# plt.legend(handles=scatter, labels=classes)
# #plt.legend((str(oneJ) for oneJ in Js_to_display), loc='upper right', fontsize=8)
# plt.savefig("LinfWRTtotalTime.eps")
# plt.savefig("LinfWRTtotalTime.jpg")
# # plt.show()

# plt.figure()
# cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
# scatter = []
# for idx,oneJ in enumerate(Js_to_display):
# 	scatter.append(plt.scatter(np.log10(computeTimePDE[oneJ]), np.log10(linferror[oneJ]), color = colours[idx], marker=markers[idx] ))
# plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
# plt.xlabel('Total PDE solve time ($log_{10}$ scale)')
# classes = ["J = " + str(oneJ) for oneJ in Js_to_display]
# plt.legend(handles=scatter, labels=classes)
# #plt.legend((str(oned) for oned in ds_to_display), loc='upper right', fontsize=8)
# plt.savefig("LinfWRTpdesolve.eps")
# plt.savefig("LinfWRTpdesolve.jpg")
# #plt.show()

# plt.figure()
# cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
# scatter = []
# for idx,oneJ in enumerate(Js_to_display):
# 	scatter.append(plt.scatter(np.log10(computeTimeRecovery[oneJ]), np.log10(linferror[oneJ]), color = colours[idx], marker=markers[idx] ))
# plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
# plt.xlabel('Total $WIHT$ recovery time ($\log_{10}$ scale)')
# classes = ["J = " + str(oneJ) for oneJ in Js_to_display]
# plt.legend(handles=scatter, labels=classes)
# #plt.legend((str(oned) for oned in ds_to_display), loc='upper right', fontsize=8)
# plt.savefig("LinfWRTl1solve.eps")
# plt.savefig("LinfWRTl1solve.jpg")
# #plt.show()

# plt.figure()
# # Display the results with J in axis
# plt.plot(Js_to_display, np.log10([computeTotalTime[j] for j in Js_to_display]), color = colours[0], marker=markers[0])
# plt.plot(Js_to_display, np.log10([computeTimePDE[j] for j in Js_to_display]), color = colours[1], marker=markers[1])
# plt.plot(Js_to_display, np.log10([computeTimeRecovery[j] for j in Js_to_display]), color = colours[2], marker=markers[2])
# plt.legend(["Total compute time", "Total PDE time", "Total sparse recovery time"])
# plt.xlabel('J (L=5)')
# plt.savefig("timesWRTjs.eps")
# plt.savefig("timesWRTjs.jpg")
# #plt.show()


fig, ax1 = plt.subplots()

ax1.set_xlabel('J')
ax1.set_ylabel('Wall clock time (h)') #, color=color)
ax1.plot(Js_to_display, [computeTotalTimeWC[d]/3600 for d in Js_to_display], color=colours[0], marker=markers[0], label='Total time')
ax1.plot(Js_to_display, [computeTotalTimeWC[d]*computeTimePDEWC[d]/3600 for d in Js_to_display], color=colours[1], marker=markers[1], label='Total PDE solve time')
ax1.plot(Js_to_display, [computeTotalTimeWC[d]*computeTimeRecoveryWC[d]/3600 for d in Js_to_display], color=colours[2], marker=markers[2], label='Total sparse recovery time')

# loc = ax1.get_xticks()
# ax1.set_xticklabels(np.arange(min(Js_to_display), max(Js_to_display) + 1, step=1))
ax1.set_xticks(range(0,max(Js_to_display)+1, 1))
# combine legends from both axes into a single stacked legend below the plot
# get handles/labels from both axes
handles1, labels1 = ax1.get_legend_handles_labels()
# ax2 may not yet have handles when called from here, so we'll get them after plotting on ax2
# ax1.tick_params(axis='y', labelcolor=color)

linf_to_display = [np.log10(linferror[d]) for d in Js_to_display]
min_linf = min(linf_to_display)
max_linf = max(linf_to_display)
delta_linf = max_linf-min_linf
ax2 = ax1.twinx() 
ax2.set_ylabel('$\log_{10}(\ell_\infty($error$))$')
ax2.plot(Js_to_display, linf_to_display, color=colours[len(Js_to_display)], marker=markers[len(Js_to_display)], label='Final approximation error')
ax2.set_ylim(min_linf-delta_linf, max_linf+delta_linf)
# ax2.set_xticks(range(0,max(Js_to_display)+1, 1))
handles2, labels2 = ax2.get_legend_handles_labels()
all_handles = handles1 + handles2
all_labels = labels1 + labels2
# place a single legend below the axes, stacked vertically (one column)
#fig.legend(all_handles, all_labels, loc='lower center', bbox_to_anchor=(0.5, -0.18), ncol=1)
fig.legend(all_handles, all_labels, loc='upper center')
# leave room at the bottom for the legend
fig.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig("timesAndErrorWRTjs.eps")
plt.savefig("timesAndErrorWRTjs.jpg")

plt.show()