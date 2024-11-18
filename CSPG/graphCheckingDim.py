from dolfin import *

#from iterative_solution import compute_true_avg_alternate as ctaa
#from iterative_solution import compute_true_avg as cta

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
import shelve

import os.path 
# thanks python >= 3.5 for the next thing! 
import pathlib # This allows to avoid the if path.exists -> mkdir thing. Hurray!

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
	else: # This is a single CSPDEResult
			pdeTimes =cspdeResultsList.t_samples
			recoveryTimes = cspdeResultsList.t_recovery
			mtxTimes = cspdeResultsList.t_matrix
			jsTimes = cspdeResultsList.t_J

	return pdeTimes, recoveryTimes

################################################3
# Meat of the script


# Keep in mind the results' folders will all have the following form: 
# Exp1H020Dim2WCosine10InfluenceTarget{StartJ}

# Remember the things? Surely there is another way to do this!
TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery', 't_J'])


# Lmax = 6
# h0 = 20
nbDim = 2
target_mesh_size = [1000]*nbDim
ds_to_display = [8, 10, 13, 16 ,20 ,25]

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



first_result	= results_all[ds_to_display[-1]] 
d 		= first_result.cspde_result[0].d # number of parameters
spde_model 	= first_result.spde_model #[TR_SPDE_MODEL_POS]
epsilon		= first_result.epsilon #[TR_EPSILON_POS]
wr_model 	= first_result.wr_model
nb_tests 	= 10 # This should be sufficient

# Get the maximum number of points in the mesh: 
#max_mesh_size = [1,1]
#for d in ds_to_display: 
#	print("Current mesh size is {} For max L = {}".format(results_all[d].spde_model.mesh_size, results_all[d].L))
#	print("All current ML results have length {}".format(len(results_all[d].cspde_result)))
#	current_target_mesh_size = 2**(len(results_all[d].cspde_result))*results_all[d].spde_model.mesh_size
#	if max_mesh_size[0] < current_target_mesh_size[0]: 
#		max_mesh_size = current_target_mesh_size

#y_GT = {}
#Z = {}
l2error = {} # np.zeros(len(ds_to_display))
linferror = {} # np.zeros(len(ds_to_display))
computeTimePDE = {} # np.zeros(len(ds_to_display))
computeTimeRecovery = {} # np.zeros(len(ds_to_display))
computeTimeMtx = {} # np.zeros(len(ds_to_display))
computeTimeJs = {} # np.zeros(len(ds_to_display))
computeTotalTime = {} # np.zeros(len(ds_to_display))
for d in ds_to_display: 
	spde_model = results_all[d].spde_model 
	spde_model.set_mesh_size(target_mesh_size)
	y_GT, Z = getGroundTruthFromModel(spde_model, wr_model, d, nb_tests, GTfolder = 'GTresults_d' + str(d) + 'H' + str(target_mesh_size[0]), GTfilenames = 'GT_')
	# y_estimated = wr_model.estimate_ML_samples(first_result[0].cspde_result, Z)
	y_estimated = results_all[idx].wr_model.estimate_ML_samples(results_all[idx].cspde_result, Z)
	l2error[d] = np.linalg.norm(y_estimated - y_GT)
	linferror[d] = np.linalg.norm(y_estimated - y_GT, ord=np.inf)
	curPdeTimes, curRecoveryTimes, curMtxTimes, curJsTimes = getComputeTimes(results_all[d].cspde_result)
	computeTimePDE[d] = curPdeTimes
	computeTimeRecovery[d] = curRecoveryTimes
	computeTimeMtx[d] = curMtxTimes
	computeTimeJs[d] = curJsTimes
	computeTotalTime[d] = curPdeTimes + curRecoveryTimes + curMtxTimes + curJsTimes

abort[0]
# Save the results: 
l2error = np.zeros(len(ds_to_display))
linferror = np.zeros(len(ds_to_display))
computeTimePDE = np.zeros(len(ds_to_display))
computeTimeRecovery = np.zeros(len(ds_to_display))
computeTimeMtx = np.zeros(len(ds_to_display))
computeTimeJs = np.zeros(len(ds_to_display))
computeTotalTime = np.zeros(len(ds_to_display))


# Let's see how our approximations perform!
for (idx, oned) in enumerate(ds_to_display): # Probably better to just read the file in this loop too instead of above.
	# Compute current estimates 
	y_estimated = results_all[idx].wr_model.estimate_ML_samples(results_all[idx].cspde_result, Z[d])
	l2error[idx] = np.linalg.norm(y_estimated - y_GT[d])
	linferror[idx] = np.linalg.norm(y_estimated - y_GT[d], ord=np.inf)
	curPdeTimes, curRecoveryTimes, curMtxTimes, curJsTimes = getComputeTimes(results_all[idx].cspde_result)
	computeTimePDE[idx] = curPdeTimes
	computeTimeRecovery[idx] = curRecoveryTimes
	computeTimeMtx[idx] = curMtxTimes
	computeTimeJs[idx] = curJsTimes
	computeTotalTime[idx] = curPdeTimes + curRecoveryTimes + curMtxTimes + curJsTimes

abort[0]
# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = np.random.randint(0, len(linferror), len(linferror)))
cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
scatter = []
for idx,oned in enumerate(ds_to_display): 
	scatter.append(plt.scatter(np.log10(computeTimePDE[idx]+computeTimeRecovery[idx]), np.log10(linferror[idx]), c = np.random.rand(3,) ))
	# scatter.append(plt.scatter(np.log10(computeTimePDE[idx]+computeTimeRecovery[idx]), np.log10(linferror[idx]), c = cmapForScatter(idx) ))
# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = [] )
plt.ylabel('$\ell_\infty$ norm of the error (via $\log_{10}$)')
plt.xlabel('Computing time ($log_{10}$ scale)')
classes = ["L = " + str(oned) for oned in ds_to_display]
plt.legend(handles=scatter, labels=classes)
#plt.legend((str(oned) for oned in ds_to_display), loc='upper right', fontsize=8)
plt.show()



