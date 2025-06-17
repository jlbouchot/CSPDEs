from dolfin import *

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
import shelve

import os.path 
# thanks python >= 3.5 for the next thing! 
import pathlib # This allows to avoid the if path.exists -> mkdir thing. Hurray!

import pandas as pd
# import seaborn as sns

from collections import namedtuple
from collections import defaultdict

from progressbar import Bar, ETA, Percentage, ProgressBar


#################################
#### This function should be externalised for reusability purposes
#################################
def getGroundTruthFromModel(spde_model, wr_model, d, nbSamples = 1000, GTfolder = 'TargetL_GTresults', GTfilenames = 'GT_'): 

	pathlib.Path(GTfolder).mkdir(parents=True, exist_ok=True) # parents = True allow to generate subfolder recursively. exist_ok prevents raising exception if the folder already exists. 
	# Folder exists for sure now!

	# Now check if a file containing sampling points exists, it is called GTfilenames +"samples.npy"
	sampleFName = GTfilenames +"samples.npy"
	if os.path.exists(os.path.join(GTfolder, sampleFName)): 
		# The file exists and we might load it!
		Z = np.load(os.path.join(GTfolder, sampleFName))
		nbExistingSamples = Z.shape[0] # [1] should be equal to the dimensionality d 
		if nbExistingSamples < nbSamples: 
			Z = np.vstack((Z,wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (nbSamples-nbExistingSamples, d)))))
	else : 
		Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (nbSamples, d) ) )
	# Now we know there are sufficiently many sampling points. We can save the file
	np.save(os.path.join(GTfolder,sampleFName), Z)

	# Step 2: generate the y values!
	y_GT = np.zeros(nbSamples)
	# Check if some of them have already been computed!
	yValuesFname = GTfilenames + "yGT.npy"
	if os.path.exists(os.path.join(GTfolder,yValuesFname)): 
		# A file exists, so let's go ahead and load it 
		existingYs = np.load(os.path.join(GTfolder,yValuesFname))
		nbExistingYs = len(existingYs)
		if nbExistingYs <= nbSamples: 
			y_GT[0:nbExistingYs] = existingYs
	else: 
		nbExistingYs = 0

	if nbExistingYs >= nbSamples: 
		return existingYs[0:nbSamples], Z[0:nbSamples]

	# Show a progressbar
	widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
	pbar    = ProgressBar(widgets=widgets)

	for k in pbar(range(nbExistingYs, nbSamples)):
		y_GT[k] = spde_model.sample(Z[k])
		# Save only every so often to avoid going to the hard memory too often
		if k+1 % 100 == 0: # k+1 because index k means we have computed k+1 data!
			np.save(os.path.join(GTfolder,yValuesFname), y_GT[0:k])

	np.save(os.path.join(GTfolder,yValuesFname), y_GT)

	return y_GT, Z


def getComputeTimes(cspdeResultsList):

	pdeTimes = 0
	recoveryTimes = 0
	mtxTimes = 0
	jsTimes = 0

	if isinstance(cspdeResultsList, list): # Run through all level 
		nbLvl = len(cspdeResultsList)
		for oneLvl in range(0,nbLvl): 
			pdeTimes = pdeTimes + cspdeResultsList[oneLvl].t_samples
			recoveryTimes = recoveryTimes + cspdeResultsList[oneLvl].t_recovery
			mtxTimes = mtxTimes + cspdeResultsList[oneLvl].t_matrix
			jsTimes = jsTimes + cspdeResultsList[oneLvl].t_J

			print(f"Level {oneLvl} has {cspdeResultsList[oneLvl].m} samples. Ansatz space has {cspdeResultsList[oneLvl].N} elements. PDE Time = {cspdeResultsList[oneLvl].t_samples}. Sparse recovery time = {cspdeResultsList[oneLvl].t_recovery}. Computing the Ansatz space took {cspdeResultsList[oneLvl].t_J}. Creating the matrix of tscheb coef took {cspdeResultsList[oneLvl].t_matrix}")
	else: # This is a single CSPDEResult
			pdeTimes =cspdeResultsList.t_samples
			recoveryTimes = cspdeResultsList.t_recovery
			mtxTimes = cspdeResultsList.t_matrix
			jsTimes = cspdeResultsList.t_J

	return pdeTimes, recoveryTimes, mtxTimes, jsTimes

################################################3
# Meat of the script


# Keep in mind the results' folders will all have the following form: 
# Exp1H020Dim2WCosine10InfluenceTarget{StartJ}

# Remember the things? Surely there is another way to do this!
TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery', 't_J'])


nbDim = 2
target_mesh_size = [3000]*nbDim
ds_to_display = [8,10,13,16, 20 ,25]

core_folder_name = 'Exp3H020Dim2WCosineDimensionalityd'
fname_to_read = 'WeightedCosine2D' # This is an unhappy mistake in my code which makes all file to have the same name. Luckily, They are all saved in separate folders. 
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
for oned in ds_to_display: 
	cur_path_to_file = core_folder_name + str(oned)
	print("Loading {0} from folder {1} ...".format(fname_to_read, cur_path_to_file))
	cur_results = sorted(shelve.open(os.path.join(cur_path_to_file,fname_to_read)).values(), key=lambda r: r.L)
	if len(cur_results) > 0: 
		# results_all.append(cur_results[0]) # This is the Check_ML.TestResult tuple
		results_all[oned] = cur_results[0]
	else: 
		# results_all.append(cur_results)
		results_all[oned] = cur_results

'''
# Debugging the load of some files: 
cur_cspde_result = results_all[8].cspde_result
for bbbb in range(len(cur_cspde_result)):
	print("Looking at J for level {} from {}".format(bbbb, len(cur_cspde_result)-1))
	print("\t J = {}".format(cur_cspde_result[bbbb].J_s))
abort[0]
'''


first_result	= results_all[ds_to_display[-1]] 
d 		= first_result.cspde_result[0].d # number of parameters
spde_model 	= first_result.spde_model #[TR_SPDE_MODEL_POS]
epsilon		= first_result.epsilon #[TR_EPSILON_POS]
wr_model 	= first_result.wr_model
nb_tests 	= 500 # This should be sufficient

# Get the maximum number of points in the mesh: 
#max_mesh_size = [1,1]
#for d in ds_to_display: 
#	print("Current mesh size is {} For max L = {}".format(results_all[d].spde_model.mesh_size, results_all[d].L))
#	print("All current ML results have length {}".format(len(results_all[d].cspde_result)))
#	current_target_mesh_size = 2**(len(results_all[d].cspde_result))*results_all[d].spde_model.mesh_size
#	if max_mesh_size[0] < current_target_mesh_size[0]: 
#		max_mesh_size = current_target_mesh_size

l2error = {} 
linferror = {} 
computeTimePDE = {} 
computeTimeRecovery = {} 
computeTimeMtx = {} 
computeTimeJs = {} 
computeTotalTime = {} 

for d in ds_to_display: 
	print("Computing results for d = {}".format(d))
	spde_model = results_all[d].spde_model 
	spde_model.set_mesh_size(target_mesh_size)
	y_GT, Z = getGroundTruthFromModel(spde_model, wr_model, 2*d, nb_tests, GTfolder = 'GTresults_d' + str(2*d) + 'H' + str(target_mesh_size[0]), GTfilenames = 'GT_')
	# y_estimated = wr_model.estimate_ML_samples(first_result[0].cspde_result, Z)
	y_estimated = results_all[d].wr_model.estimate_ML_samples(results_all[d].cspde_result, Z)
	l2error[d] = np.linalg.norm(y_estimated - y_GT)
	linferror[d] = np.linalg.norm(y_estimated - y_GT, ord=np.inf)
	curPdeTimes, curRecoveryTimes, curMtxTimes, curJsTimes = getComputeTimes(results_all[d].cspde_result)
	computeTotalTime[d] = curPdeTimes + curRecoveryTimes + curMtxTimes + curJsTimes
	computeTimePDE[d] = curPdeTimes / computeTotalTime[d]
	computeTimeRecovery[d] = curRecoveryTimes / computeTotalTime[d]
	computeTimeMtx[d] = curMtxTimes / computeTotalTime[d]
	computeTimeJs[d] = curJsTimes / computeTotalTime[d]
	
all_data_df = pd.DataFrame(data=[computeTimePDE, computeTimeRecovery, computeTimeMtx, computeTimeJs, computeTotalTime]).rename({0: 'PDE', 1: 'Recovery', 2: 'Matrix', 3: 'J', 4: 'Total'}).transpose()

colours = plt.cm.rainbow(np.linspace(0, 1, len(ds_to_display)))
markers = ["o", "v", "^", "<", ">", "s", "8"]


# # We have here the graph of compute time wrt target accuracy, for various dimensions 
# cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
# plt.figure()
# scatter = []
# for idx, oned in enumerate(ds_to_display): 
# 	scatter.append(plt.scatter(np.log10(computeTotalTime[oned]), np.log10(linferror[oned]), color = colours[idx], marker=markers[idx] ))
# plt.ylabel('$\log_{10}(\ell_\infty($error$))$')
# plt.xlabel('Computing time ($log_{10}$ scale)')
# classes = ["d = " + str(oned) for oned in ds_to_display]
# plt.legend(handles=scatter, labels=classes)
# #plt.legend((str(oned) for oned in ds_to_display), loc='upper right', fontsize=8)
# plt.show()

# # # Here we have the ratio of computeTime vs ratio of recovery time
# # cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
# # scatter = []
# # plt.figure()
# # for idx, oned in enumerate(ds_to_display): 
# # 	scatter.append(plt.scatter(computeTimePDE[oned], computeTimeRecovery[oned], color = colours[idx], marker=markers[idx] ))
# # plt.ylabel('Fraction of time used for sparse recovery')
# # plt.xlabel('Fraction of time used for PDE computation')
# # classes = ["d = " + str(oned) for oned in ds_to_display]
# # plt.legend(handles=scatter, labels=classes)
# # #plt.legend((str(oned) for oned in ds_to_display), loc='upper right', fontsize=8)
# # plt.show()

# # We may now disply as bar graphs the ratios of various times with changing dimensions of the data
# plt.figure()
# plt.bar(x = ds_to_display, height=np.array(all_data_df[["PDE"]]).transpose()[0])
# plt.bar(x = ds_to_display, height=np.array(all_data_df[["Recovery"]]).transpose()[0], bottom=np.array(all_data_df[["PDE"]]).transpose()[0])
# plt.bar(x = ds_to_display, height=np.array(all_data_df[["Matrix"]]).transpose()[0], bottom=np.array(all_data_df[["PDE"]]).transpose()[0]+np.array(all_data_df[["Recovery"]]).transpose()[0])
# plt.bar(x = ds_to_display, height=np.array(all_data_df[["J"]]).transpose()[0], bottom=np.array(all_data_df[["PDE"]]).transpose()[0]+np.array(all_data_df[["Recovery"]]).transpose()[0]+np.array(all_data_df[["J"]]).transpose()[0])
# plt.legend(["PDE","Sparse recovery","Matrix computation","Index set construction"])
# plt.show()

# # 4 Graphs representing the four curves with respect to the changing dimensionality
# plt.figure()
# plt.plot(all_data_df[["PDE","Recovery","Matrix","J"]])
# plt.legend(["PDE","Sparse recovery","Matrix computation","Index set construction"])
# plt.show()



# 5 All on the same graph!
colours = plt.cm.rainbow(np.linspace(0, 1, len(ds_to_display)))
markers = ["o", "v", "^", "<", ">", "s", "8"]

ds = [2*d for d in ds_to_display]

fig, ax1 = plt.subplots()

ax1.set_xlabel('d')
ax1.set_ylabel('time (h)') #, color=color)
ax1.plot(ds, [computeTotalTime[d]/3600 for d in ds_to_display], color=colours[0], label='Total time')
ax1.plot(ds, [computeTotalTime[d]*computeTimePDE[d]/3600 for d in ds_to_display], color=colours[1], label='Total PDE solve time')
ax1.plot(ds, [computeTotalTime[d]*computeTimeRecovery[d]/3600 for d in ds_to_display], color=colours[2], label='Total sparse recovery time')
ax1.set_xticklabels(ax1.get_xticks().astype(int))
ax1.legend()
# ax1.tick_params(axis='y', labelcolor=color)

linf_to_display = [np.log10(linferror[d]) for d in ds_to_display]
min_linf = min(linf_to_display)
max_linf = max(linf_to_display)
delta_linf = max_linf-min_linf
ax2 = ax1.twinx() 
ax2.set_ylabel('$\log_{10}(\ell_\infty($error$))$')
ax2.plot(ds, linf_to_display, color=colours[3], label='Final approximation error')
ax2.set_ylim(min_linf-delta_linf, max_linf+delta_linf)
ax2.legend()
fig.tight_layout()

plt.savefig("timesAndErrorWRTDs.eps")
plt.savefig("timesAndErrorWRTDs.jpg")

plt.show()