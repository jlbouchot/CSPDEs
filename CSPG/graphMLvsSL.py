from dolfin import *

from iterative_solution import compute_true_avg_alternate as ctaa
from iterative_solution import compute_true_avg as cta

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

from progressbar import Bar, ETA, Percentage, ProgressBar


#################################
#### This function should be externalised for reusability purposes
#################################
def getGroundTruthFromModel(spde_model, wr_model, d, nbSamples = 1000, GTfolder = 'GTresults', GTfilenames = 'GT_'): 

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

	if isinstance(cspdeResultsList, list): # Run through all level 
		nbLvl = len(cspdeResultsList)
		for oneLvl in range(0,nbLvl): 
			pdeTimes = pdeTimes + cspdeResultsList[oneLvl].t_samples
			recoveryTimes = recoveryTimes + cspdeResultsList[oneLvl].t_recovery
	else: # This is a single CSPDEResult
			pdeTimes =cspdeResultsList.t_samples
			recoveryTimes = cspdeResultsList.t_recovery

	return pdeTimes, recoveryTimes

################################################3
# Meat of the script


# Keep in mind the results' folders will all have the following form: 
# Exp1Dim2WCosine20InfluenceJ0

# Remember the things? Surely there is another way to do this!
TestResult = namedtuple('TestResult', ['spde_model', 'wr_model', 'epsilon', 'L', 'cspde_result'])
CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery'])


Lmax = 6
h0 = 10;
Js_to_display = [1,2,3,4,5,6] # Note that J = 6 corresponds to the SL appraoch

core_folder_name = 'Exp1Dim2WCosine20InfluenceJ'
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

# Load and plot results for all J's, one after the other 
results_all = []
for oneJ in Js_to_display: 
	cur_path_to_file = core_folder_name + str(oneJ)
	print("Loading {0} from folder {1} ...".format(fname_to_read, cur_path_to_file))
	cur_results = sorted(shelve.open(os.path.join(cur_path_to_file,fname_to_read)).values(), key=lambda r: r.L)
	if len(cur_results) > 0: 
		results_all.append(cur_results[0]) # This is the Check_ML.TestResult tuple
	else: 
		results_all.append(cur_results)



first_result	= results_all[:-1] # At that moment, first_result is a single-level result
d 		= first_result[0].cspde_result[0].d # number of parameters
#d            = first_result[0][TR_CR_POS][CR_D_POS] # number of parameters
spde_model 	= first_result[0].spde_model #[TR_SPDE_MODEL_POS]
epsilon		= first_result[0].epsilon #[TR_EPSILON_POS]
wr_model 	= first_result[0].wr_model
nb_tests 	= 100#00 # This should be sufficient


y_GT, Z = getGroundTruthFromModel(spde_model, wr_model, d, nb_tests, GTfolder = 'GTresults', GTfilenames = 'GT_')
y_estimated = wr_model.estimate_ML_samples(first_result[0].cspde_result, Z)


# Save the results: 
l2error = np.zeros(len(Js_to_display))
linferror = np.zeros(len(Js_to_display))
computeTimePDE = np.zeros(len(Js_to_display))
computeTimeRecovery = np.zeros(len(Js_to_display))


# Let's see how our approximations perform!
for (idx, oneJ) in enumerate(Js_to_display): # Probably better to just read the file in this loop too instead of above.
	# Compute current estimates 
	y_estimated = results_all[idx].wr_model.estimate_ML_samples(results_all[idx].cspde_result, Z)
	l2error[idx] = np.linalg.norm(y_estimated - y_GT)
	linferror[idx] = np.linalg.norm(y_estimated - y_GT, ord=np.inf)
	curPdeTimes, curRecoveryTimes = getComputeTimes(results_all[idx].cspde_result)
	computeTimePDE[idx] = curPdeTimes
	computeTimeRecovery[idx] = curRecoveryTimes


# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = np.random.randint(0, len(linferror), len(linferror)))
cmapForScatter = plt.cm.get_cmap('hsv', len(linferror))
scatter = []
for idx,oneJ in enumerate(Js_to_display): 
	scatter.append(plt.scatter(np.log10(computeTimePDE[idx]+computeTimeRecovery[idx]), np.log10(linferror[idx]), c = np.random.rand(3,) ))
	# scatter.append(plt.scatter(np.log10(computeTimePDE[idx]+computeTimeRecovery[idx]), np.log10(linferror[idx]), c = cmapForScatter(idx) ))
# scatter=plt.scatter(np.log10(computeTimePDE+computeTimeRecovery), np.log10(linferror), c = [] )
plt.ylabel('$\ell_\infty$ norm of the error (via $\log_10$)')
plt.xlabel('Computing time ($log_10$ scale)')
classes = ["J = " + str(oneJ) for oneJ in Js_to_display]
plt.legend(handles=scatter, labels=classes)
#plt.legend((str(oneJ) for oneJ in Js_to_display), loc='upper right', fontsize=8)
plt.show()



