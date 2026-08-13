# A utility package for all sorts of helper functions for MLCSPG

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019 - 2026, INRIA, LMU Munisch, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__create__ = "2026/01/22"
__lastmodified__ = "2026/08/07"

import WR
import resource as res
import time
import numpy as np
import pathlib
import os.path
from progressbar import Bar, ETA, Percentage, ProgressBar

def get_sampling_type(sampling_name):
    switcher = {
        "pragmatic": WR.cs_pragmatic_m,
        "theoretic": WR.cs_theoretic_m,
        "new": WR.cs_theoretic_m_new,
        "p": WR.cs_pragmatic_m,
        "t": WR.cs_theoretic_m,
    }
    return switcher.get(sampling_name, WR.cs_pragmatic_m)

def is_tensor_based_possible(method):
    """ 
    Returns True if the method is compatible with tensor-based operators, False otherwise. 
    It is false for anything based on convex optimization, as these methods are not compatible with tensor-based operators.
    """
    return method not in ["bp", "bpdn"]

def time_things(start_time = None):
    """
    Returns the time elapsed since start_time, in seconds. If start_time is None, returns the current time and resource usage.
    Note that it returns the wall-clock time, user CPU time and system CPU time.
    Total CPU time may be obtained by adding user and system CPU times.
    """
    tWC = time.time()
    start = res.getrusage(res.RUSAGE_SELF)
    if start_time: 
        tWC = tWC - start_time[0]
        tUser = start.ru_utime - start_time[1].ru_utime
        tSys = start.ru_stime - start_time[1].ru_stime
        return tWC, tUser, tSys
    else:
        return tWC, start
    # tCPU = tUser + tSys

def getComputeTimes(cspdeResultsList):

    pdeTimes = [0]*3 # Contains WC, user, system times
    recoveryTimes = [0]*3
    mtxTimes = [0]*3
    jsTimes = [0]*3

    if isinstance(cspdeResultsList, list): # Run through all level 
        nbLvl = len(cspdeResultsList)
        for oneLvl in range(0,nbLvl):
            pdeTimes = [x + y for x, y in zip(pdeTimes, cspdeResultsList[oneLvl].t_samples)]
            recoveryTimes = [x + y for x, y in zip(recoveryTimes, cspdeResultsList[oneLvl].t_recovery)]
            mtxTimes = [x + y for x, y in zip(mtxTimes, cspdeResultsList[oneLvl].t_matrix)]
            # WARNING: This is a bug corrected for the long run experiments. Right now, we need to duplicate the results in some cases.
            # REMOVE: print(f"Type of t_J is {type(cspdeResultsList[oneLvl].t_J)}")
            if (type(cspdeResultsList[oneLvl].t_J) is list) or (type(cspdeResultsList[oneLvl].t_J) is np.ndarray) or (type(cspdeResultsList[oneLvl].t_J) is tuple):
                # print(f"t_J is a list of length {len(cspdeResultsList[oneLvl].t_J)}")
                jsTimes = [x + y for x, y in zip(jsTimes, cspdeResultsList[oneLvl].t_J)]
            else:
                # print(f"t_J is a single value {cspdeResultsList[oneLvl].t_J}")

                jsTimes = [x + y for x, y in zip(jsTimes, [cspdeResultsList[oneLvl].t_J]*3)]
    else: # This is a single CSPDEResult
            pdeTimes =cspdeResultsList.t_samples
            recoveryTimes = cspdeResultsList.t_recovery
            mtxTimes = cspdeResultsList.t_matrix
            if type(cspdeResultsList[oneLvl].t_J) is list:
                jsTimes = cspdeResultsList.t_J
            else:
                jsTimes = [cspdeResultsList.t_J]*3

    return pdeTimes, recoveryTimes, mtxTimes, jsTimes

def getGroundTruthFromModel(spde_model, wr_model, d, nbSamples = 1000, GTfolder = 'TargetL_GTresults', GTfilenames = 'GT_'): 

    pathlib.Path(GTfolder).mkdir(parents=True, exist_ok=True) # parents = True allow to generate subfolder recursively. exist_ok prevents raising exception if the folder already exists. 
    # Folder exists for sure now!

    # Now check if a file containing sampling points exists, it is called GTfilenames +"samples.npy"
    sampleFName = GTfilenames +"samples.npy"
    print(f"Checking if {os.path.join(GTfolder, sampleFName)} exists before relaunching computations...")
    if os.path.exists(os.path.join(GTfolder, sampleFName)): 
        print(f"Found {os.path.join(GTfolder, sampleFName)}. Loading it...")
        # The file exists and we might load it!
        Z = np.load(os.path.join(GTfolder, sampleFName))
        nbExistingSamples = Z.shape[0] # [1] should be equal to the dimensionality d 
        if nbExistingSamples < nbSamples: 
            Z = np.vstack((Z,wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (nbSamples-nbExistingSamples, d)))))
    else : 
        print(f"File {os.path.join(GTfolder, sampleFName)} does not exist. Computing samples...")
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

    print(f"Found {nbExistingYs} existing samples out of the expected {nbSamples}")
    if nbExistingYs >= nbSamples: 
        return existingYs[0:nbSamples], Z[0:nbSamples]

    # Show a progressbar
    # spde_model.refine_mesh(2)
    print("Current FEM width is {}".format(spde_model.mesh_size[0]))
    widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
    pbar    = ProgressBar(widgets=widgets)
    
    gt_params = {"linear_solver": "gmres", "preconditioner": "amg", "relative_tolerance": 1e-7, "absolute_tolerance": 1e-12}

    for k in pbar(range(nbExistingYs, nbSamples)):
        y_GT[k] = spde_model.sample(Z[k], gt_params)
        # Save only every so often to avoid going to the hard memory too often
        if ((k+1) % 100) == 0: # k+1 because index k means we have computed k+1 data!
            np.save(os.path.join(GTfolder,yValuesFname), y_GT[0:k])

    np.save(os.path.join(GTfolder,yValuesFname), y_GT)

    return y_GT, Z