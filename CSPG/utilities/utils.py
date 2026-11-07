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
__lastmodified__ = "2026/11/07"

import WR
import SPDESampler.Functionals as F
import numpy as np

from utilities import gt_meta
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

def get_operator_type(operator_name):
    switcher = {
        "Chebyshev": WR.Chebyshev,
        "Cheb_Alt": WR.Cheb_Alt,
        "Cheb_Sparse": WR.Cheb_Sparse,
        "Cheb_AltPP": WR.Cheb_AltPP,
    }
    return switcher.get(operator_name, WR.Cheb_Sparse)

def get_functional_type(functional_name):
    '''
    Select the appropriate functional based on the experiment in progress.
    WARNING: This function needs to be updated if new functionals are added to the SPDE module. It is a bit of a hack, but it works for now.
    TODO: Eventually, we should find a better way to handle this, maybe by using a dictionary or a class that maps functional names to their corresponding classes.
    TODO: Make sure that parametric functionals are handled correctly. 
    TODO: Question for future: Could the functional be part of the parametric PDE too? 
    '''
    switcher = {
        "average": F.AverageQOI,
        "avg": F.AverageQOI,
        "box": F.BoxAverage,
    }
    ## The entries are CLASSES, not instances: a dict literal of instances would
    ## build every functional on each call, and an unknown name must fail loudly
    ## rather than silently fall back to a different quantity of interest.
    if functional_name not in switcher:
        raise ValueError("Unknown quantity of interest '{0}'. Available: {1}".format(
            functional_name, sorted(switcher)))
    return switcher[functional_name]()

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

def loadGroundTruth(d, nbSamples = 1000, GTidentity = None, GTfolder = None, GTfilenames = 'GT_'):
    """
    Load the ground truth that belongs to GTidentity, or say exactly how to compute it.

    If too few samples are available it tells how to compute more.

    Generating the reference is the job of one of

        testing/compute_ground_truth.py     (one mesh)
        testing/richardson_extrapolation.py (extrapolated in h)

    which record their solver and their convergence in GT_meta.json. This
    function only reads, and refuses to read a set that was computed for a
    different problem. Build GTidentity with gt_meta.identity_from_cfg.
    """
    if GTidentity is None:
        raise ValueError(
            "loadGroundTruth needs GTidentity: the folder name alone does not "
            "say which coefficient, QoI or seed the reference belongs to. Use "
            "gt_meta.identity_from_cfg(<experiment .ini>, mesh, seed) and set its 'd' if the "
            "experiment sweeps d.")

    folder = GTfolder or gt_meta.folder(GTidentity)
    sampleFName = os.path.join(folder, GTfilenames + "samples.npy")
    yValuesFname = os.path.join(folder, GTfilenames + "yGT.npy")

    ## One place says what to run, so the message is the same wherever the
    ## reference turns out to be missing or short.
    howto = ("Run the ground truth computations with the same .ini, e.g.\n"
             "    python compute_ground_truth.py --cfg <experiment .ini> "
             "--mesh {0} --nb_samples {1} --seed {2}".format(
                 GTidentity["mesh"], nbSamples, GTidentity["seed"]))

    if not os.path.exists(yValuesFname):
        raise SystemExit("no ground truth in {0}\n  identity: {1}\n{2}".format(
            folder, gt_meta.canonical(GTidentity), howto))
    gt_meta.check(folder, GTidentity)       # exits if that folder holds another problem

    ## The file is nb_samples long with NaN where a sample is still missing,
    ## so completeness is checked, not inferred from the length.
    y_GT = np.load(yValuesFname)
    if len(y_GT) < nbSamples:
        raise SystemExit(
            "{0} holds {1} values but {2} were asked for.\n{3}".format(
                yValuesFname, len(y_GT), nbSamples, howto))
    hole = np.flatnonzero(~np.isfinite(y_GT[:nbSamples]))
    if hole.size:
        raise SystemExit(
            "{0} is missing {1} of the first {2} values (samples {3}{4}).\n{5}".format(
                yValuesFname, hole.size, nbSamples, hole[:20].tolist(),
                " ..." if hole.size > 20 else "", howto))

    ## The points are regenerated rather than trusted: gt_meta.draw is the one
    ## definition of the sampling measure, and the stored file is only a
    ## cross-check that the reference really used it.
    Z = gt_meta.draw(nbSamples, d, GTidentity["seed"])
    if os.path.exists(sampleFName):
        stored = np.load(sampleFName)
        k = min(stored.shape[0], nbSamples)
        if stored.shape[1] != d or not np.allclose(stored[:k], Z[:k]):
            raise SystemExit(
                "{0} does not hold the points gt_meta.draw(seed={1}) produces: the "
                "reference was generated with a different measure or a different d. "
                "Move it aside and regenerate.".format(sampleFName, GTidentity["seed"]))

    print("ground truth: {0} values from {1}".format(nbSamples, folder))
    meta = gt_meta.read(folder) or {}
    q = meta.get("quality", {})
    if q:
        print("  solver {0}, rtol {1}, all_converged {2}".format(
            q.get("solver"), q.get("rtol"), q.get("all_converged")))
        if q.get("all_converged") is False:
            print("  WARNING: that set has non-converged samples: {0}".format(
                q.get("failed_samples")))
    return y_GT[0:nbSamples], Z[0:nbSamples]


def getGroundTruthFromModel(spde_model, wr_model, d, nbSamples = 1000, GTfolder = None, GTfilenames = 'GT_', GTidentity = None):
    """Deprecated: the model arguments are no longer used, because nothing is
    computed here any more. Kept so an older script still runs; call
    loadGroundTruth directly."""
    return loadGroundTruth(d, nbSamples, GTidentity, GTfolder, GTfilenames)

def default_tolerance_sparse_recovery(targetL, nb_recoveries, h0, t=1, tprime=1, constant=0.125):
    """
    Returns the default tolerance for sparse recovery algorithms. 
    This is a function to avoid hardcoding the value in multiple places.
    #TODO: 
    Add explanation about the parameters and the formula. This is a function of the target level, the number of recoveries, the initial mesh size, and some parameters t, tprime and constant.

    """
    return 1/(h0*h0*4**targetL)/(2**(t+tprime)-1)/nb_recoveries*constant