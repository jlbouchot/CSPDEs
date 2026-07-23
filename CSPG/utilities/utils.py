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
__lastmodified__ = "2026/07/22"

import WR
import resource as res
import time

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
