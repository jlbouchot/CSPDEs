from dolfin import *

#from iterative_solution import compute_true_avg_alternate as ctaa
#from iterative_solution import compute_true_avg as cta

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
__lastmodified__ = "2026/08/07"

# Keep in mind the results' folders will all have the following form: 
# TODO: Add proper naming conventions

# Remember the things? Surely there is another way to do this! ==> Pass it in the utils? 
TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery'])


Lmax = [2,3,4]
nbDim = 2
target_mesh_size = [3000]*nbDim # TODO: Will need to change this!!!
Js_to_display = [l - 2 for l in Lmax] # Note that J = 6 corresponds to the SL appraoch


core_folder_name = os.path.join('results', 'Exp1_debug')
base_fname = 'L_val_'
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

# Load and plot results for all J's, one after the other 
results_all = {}
for idx, oneJ in enumerate(Js_to_display): 
	cur_path_to_file = core_folder_name # + str(oneJ)
	fname_to_read = base_fname + str(Lmax[idx])
	print("Loading {0} from folder {1} ...".format(fname_to_read, cur_path_to_file))
	cur_results = sorted(shelve.open(os.path.join(cur_path_to_file,fname_to_read)).values(), key=lambda r: r.L)
	# print(cur_results)
	if len(cur_results) > 0: 
		results_all[oneJ] = cur_results[0] # This is the Check_ML.TestResult tuple
	else: 
		results_all[oneJ] = cur_results
	# print("Current results has length {} while its first element has a mesh size of {}".format(len(results_all[:-1]), results_all[:-1][0].spde_model.mesh_size[0] ))
	# print("Current results has length {}".format(len(results_all) ))
	#print("Results all at 0 is {}".format(results_all[0]))
	#print("Results all is {}".format(results_all))
	#print("Last element in results all is ".format(results_all[:-1]))
	#print("Type of last element of results all is {}".format(type(results_all[:-1])))
	# print("Length of the last element of results all is {}".format(len(results_all[:-1])))
	# print(results_all[-1][0])
	# print("Final discretization mesh is{}".format(results_all[-1].spde_model.mesh_size[0]))
	# print("Latest results is list of length {}".format(len(results_all[-1]) ))
	# print("Nb of levels should be 3 == {} ?".format(len(results_all[-1].cspde_result)))



finest_result	= results_all[Js_to_display[-1]]
d 		= finest_result.cspde_result[0].d # number of parameters
#d            = finest_result[0][TR_CR_POS][CR_D_POS] # number of parameters
spde_model 	= finest_result.spde_model #[TR_SPDE_MODEL_POS]
# epsilon		= finest_result.epsilon #[TR_EPSILON_POS]
wr_model 	= finest_result.wr_model
nb_tests 	= 1000

spde_model.set_mesh_size(target_mesh_size)
print("Getting ground truth from model with width {}".format(spde_model.mesh_size[0]))
y_GT, Z = getGroundTruthFromModel(spde_model, wr_model, d, nb_tests, GTfolder = 'GTresults_d' + str(d) + 'H' + str(target_mesh_size[0]), GTfilenames = 'GT_')
# y_estimated = wr_model.estimate_ML_samples(finest_result[0].cspde_result, Z)


# # Save the results: 
# l2error = np.zeros(len(Js_to_display))
# linferror = np.zeros(len(Js_to_display))
# computeTimePDE = np.zeros(len(Js_to_display))
# computeTimeRecovery = np.zeros(len(Js_to_display))

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

# Let's see how our approximations perform!
for oneJ in Js_to_display: # Probably better to just read the file in this loop too instead of above.
	# Compute current estimates 
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

	# l2error[idx] = np.linalg.norm(y_estimated - y_GT)
	# linferror[idx] = np.linalg.norm(y_estimated - y_GT, ord=np.inf)
	# curPdeTimes, curRecoveryTimes = getComputeTimes(results_all[idx].cspde_result)
	# computeTimePDE[idx] = curPdeTimes
	# computeTimeRecovery[idx] = curRecoveryTimes

all_data_df = pd.DataFrame(data=[computeTimePDEWC, computeTimeRecoveryWC, computeTimeMtxWC, computeTimeJsWC, computeTotalTimeWC, computeTimePDEUser, computeTimeRecoveryUser, computeTimeMtxUser, computeTimeJsUser, computeTotalTimeUser, computeTimePDESys, computeTimeRecoverySys, computeTimeMtxSys, computeTimeJsSys, computeTotalTimeSys]).rename({0: 'PDE_WC', 1: 'Recovery_WC', 2: 'Matrix_WC', 3: 'J_WC', 4: 'Total_WC', 5: 'PDE_User', 6: 'Recovery_User', 7: 'Matrix_User', 8: 'J_User', 9: 'Total_User', 10: 'PDE_Sys', 11: 'Recovery_Sys', 12: 'Matrix_Sys', 13: 'J_Sys', 14: 'Total_Sys'}).transpose()


colours = plt.cm.rainbow(np.linspace(0, 1, len(Js_to_display)))
markers = ["o", "v", "^", "<", ">", "s", "8"]

# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = np.random.randint(0, len(linferror), len(linferror)))
plt.figure()
cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
scatter = []
for idx, oneJ in enumerate(Js_to_display): 
	scatter.append(plt.scatter(np.log10(computeTotalTimeWC[oneJ]), np.log10(linferror[oneJ]), c = colours[idx], marker=markers[idx] ))
	#scatter.append(plt.scatter(np.log10(computeTimePDEWC[oneJ]+computeTimeRecoveryWC[oneJ]), np.log10(linferror[oneJ]), c = np.random.rand(3,) ))
	# scatter.append(plt.scatter(np.log10(computeTimePDEWC[idx]+computeTimeRecoveryWC[idx]), np.log10(linferror[idx]), c = cmapForScatter(idx) ))
# scatter=plt.scatter(np.log10(computeTimePDEWC+computeTimeRecoveryWC), np.log10(linferror), c = [] )
plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
plt.xlabel('Wall Clock Computing time ($log_{10}$ scale)')
classes = ["L = " + str(oneJ) for oneJ in Js_to_display]
plt.legend(handles=scatter, labels=classes)
#plt.legend((str(oneJ) for oneJ in Js_to_display), loc='upper right', fontsize=8)
plt.show()

plt.figure()
cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
scatter = []
for idx, oneJ in enumerate(Js_to_display): 
	scatter.append(plt.scatter(2**(-oneJ)*target_mesh_size[0], linferror[oneJ], c = colours[idx], marker=markers[idx] ))
	# scatter.append(plt.scatter(2**(-oneJ)*target_mesh_size[0], np.log10(linferror[oneJ]), c = np.random.rand(3,) ))
	#scatter.append(plt.scatter(np.log10(computeTimePDE[oneJ]+computeTimeRecovery[oneJ]), np.log10(linferror[oneJ]), c = np.random.rand(3,) ))
	# scatter.append(plt.scatter(np.log10(computeTimePDE[idx]+computeTimeRecovery[idx]), np.log10(linferror[idx]), c = cmapForScatter(idx) ))
# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = [] )
plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
plt.xlabel('Target accuracy ($2^{-L}h_0$)')
classes = ["L = " + str(oneJ) for oneJ in Js_to_display]
plt.legend(handles=scatter, labels=classes)
#plt.legend((str(oneJ) for oneJ in Js_to_display), loc='upper right', fontsize=8)
plt.show()



fig, ax1 = plt.subplots()

ax1.set_xlabel('L')
ax1.set_ylabel('time (h)') #, color=color)
ax1.plot(Js_to_display, [computeTotalTimeWC[d]/3600 for d in Js_to_display], color=colours[0], marker=markers[0], label='Total time')
ax1.plot(Js_to_display, [computeTotalTimeWC[d]*computeTimePDEWC[d]/3600 for d in Js_to_display], color=colours[1], marker=markers[1], label='Total Wall Clock PDE solve time (hours)')
ax1.plot(Js_to_display, [computeTotalTimeWC[d]*computeTimeRecoveryWC[d]/3600 for d in Js_to_display], color=colours[2], marker=markers[2], label='Total Wall Clock sparse recovery time (hours)')
ax1.set_xticks(range(0,max(Js_to_display)+1, 1))
ax1.legend()
# ax1.tick_params(axis='y', labelcolor=color)

linf_to_display = [np.log10(linferror[d]) for d in Js_to_display]
min_linf = min(linf_to_display)
max_linf = max(linf_to_display)
delta_linf = max_linf-min_linf
ax2 = ax1.twinx() 
ax2.set_ylabel('$\log_{10}(\ell_\infty($error$))$')
ax2.plot(Js_to_display, linf_to_display, color=colours[3], marker=markers[3], label='Final approximation error')
ax2.set_ylim(min_linf-delta_linf, max_linf+delta_linf)
# ax2.set_xticks(range(0,max(Js_to_display)+1, 1))
ax2.legend()
fig.tight_layout()
plt.savefig("timesAndErrorWRTLs.eps")
plt.savefig("timesAndErrorWRTLs.jpg")

plt.show()
