import numpy as np
import itertools
from collections import namedtuple

import WR.Operators.Chebyshev as Chebyshev
import WR.Operators.operator_from_matrix as ofm

import time

import os.path

## Add utilities to the path
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'utilities'))
import utils as u

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019 - 2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/10/09"

from time import sleep

CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery', 't_J'])

def level_plan(cfg):
    """Return [(level, s_level, is_first_level)] exactly as in CSPDE_ML."""
    L_first, L = cfg["l_start"], cfg["nb_level"]
    p, p0 = cfg["p_t"], cfg["p_0"]
    rate = cfg["t_0"] + cfg["t_prime"]

    s_L = np.ceil((cfg["dat_constant"] * (L - L_first)) ** (p / (1 - p)))
    s_J = np.ceil(cfg["const_sj"] ** (p0 / (1 - p0)) * 2 ** (L * p0 * rate / (1 - p0)))
    s_J = max(np.ceil(s_L * 2 ** ((L - L_first) * rate * p / (1 - p))), s_J)

    plan = [(L_first, s_J, True)]
    plan += [(l, np.ceil(s_L * 2 ** ((L - l) * rate * p / (1 - p))), False) for l in range(L_first + 1, L + 1)]
    return plan

def _plan_level(level, s, is_first, ctx):
    """Index set and sample size of a level: cheap, done for all levels before any PDE solve."""
    wr_model = ctx["wr"]
    ctx["log"]("Generating an Ansatz space of multiindices " +  ("that have total degree <= {}".format(max_degree) if ctx["ansatz_space"] else "of weigthed sparsity {}".format(s) ) )
    J_s, t_J = index_set(s, wr_model, ctx["ansatz_space"]) ## TODO: Check function call
    N, d = len(J_s), len(J_s[0])
    m = wr_model.get_m_from_s_N(s, N)
    ctx["log"]("Level {0} ({1}): s = {2}, N = {3}, m = {4}, d = {5}, epsilon = {6}".format(
        level, "single level approximation" if is_first else "detail", s, N, m, d, ctx["tol_res"] * np.sqrt(m)))
    wr_model.check(N, m)
    return dict(level=level, s=s, is_first=is_first, J_s=J_s, t_J=t_J, N=N, d=d, m=m)


def _validate_plan(plan):
    # s <= 4 gives J_s = {0}: N = 1 and m = ceil(2 s log 1) = 0. The legacy code then divides by
    # sqrt(m) = 0 in the operator and crashes, possibly after hours spent on the coarser levels.
    bad = [p for p in plan if p["m"] < 1]
    if bad:
        raise ValueError("No samples would be drawn at level(s) {0} (s = {1}, N = {2}). Increase dat_constant / const_sj "
                         "or reduce nb_level - l_start.".format([p["level"] for p in bad], [p["s"] for p in bad], [p["N"] for p in bad]))


def CSPDE_ML(spde_model, wr_model, dict_config, sparse_config): 
    """Multi-level compressed sensing solver for SPDE approximation.

    This routine coordinates a multi-level compressed-sensing workflow for
    approximating a parametric SPDE. For each level determined by
    :func:`level_plan`, it builds an index set, draws samples of the PDE
    solution, assembles (or loads) the sampling matrix, computes recovery
    weights, and performs a weighted sparse recovery using the provided
    weighted recovery model (`wr_model`). Results for each level are
    returned as a list of ``CSPDEResult`` namedtuples.

    Parameters
    ----------
    spde_model : object
        SPDE model object exposing at least the methods ``samples(Z)`` and
        ``refine_mesh()``. ``samples(Z)`` should return a vector of
        evaluations of the quantity of interest at the sample points ``Z``.
    wr_model : object
        Weighted-recovery model that provides the sampling operator
        interface through ``wr_model.operator`` and a recovery method
        ``wr_model.method(A, y, w, s, epsilon, nb_iter, print_every)``. It
        must also provide ``get_m_from_s_N(s, N)`` and ``check(N, m)``.
    dict_config : dict
        Global configuration dictionary. Expected keys include
        ``nb_level``, ``l_start``, ``ansatz_space``, ``experiment_name``,
        ``do_tensor`` and optional ``seed`` and ``no_compute``.
    sparse_config : dict
        Sparse-recovery configuration. Expected keys include ``tol_res``,
        ``nb_iter`` and optional ``log`` (logging frequency / verbosity).

    Returns
    -------
    list of CSPDEResult
        A list containing one ``CSPDEResult`` namedtuple per processed
        level. Each entry contains fields
        (``J_s, N, s, m, d, Z, y, A, w, result, t_samples, t_matrix,
        t_recovery, t_J``) describing the polynomial index set, problem
        sizes, sampling points, measurements, sampling operator (or 0 if
        not stored), weights, recovery result object, and timing
        information.

    Notes
    -----
    - The function may read/write cached sampling points and sampling
      matrices under the directory given by ``dict_config['experiment_name']``.
    - If ``dict_config['no_compute']`` is True, an empty list is returned.
    - Randomness is controlled by ``dict_config['seed']`` when present.
    - Side effects: prints progress messages and may call
      ``spde_model.refine_mesh()`` when advancing levels.
    """
    # First set up some basic things needed for the rest of the computations
    lvl_by_lvl_result = [] # This will keep the results

    log_freq = sparse_config.get("log", None) # TODO: Eventually add other logging tools, in which case this parameter name will no longer be acceptable
    seed = dict_config.get("seed", None)
    ctx = dict(
        spde=spde_model, wr=wr_model, rng=np.random.default_rng(None if seed is None else int(seed)), seed=seed,
        L=dict_config["nb_level"], L_first=dict_config["l_start"], ansatz_space=dict_config["ansatz_space"],
        no_compute=dict_config["no_compute"], tol_res=sparse_config["tol_res"], nb_iter=sparse_config["nb_iter"],
        log_freq=log_freq, log=(lambda msg: print(msg)),
    )

    sampling_fname = os.path.join(dict_config["experiment_name"], 'sampling_points_Lmax' + str(ctx["L"]))
    # datamtx_fname = os.path.join(dict_config["experiment_name"], 'datamtx' + ('_tensor' if dict_config["do_tensor"] else '_no_tensor') + '_Lmax' + str(ctx["L"]))
    datamtx_fname = os.path.join(dict_config["experiment_name"], 'datamtx' + '_' +dict_config["operator"] + '_Lmax' + str(ctx["L"]))

    plan = [_plan_level(level, s, is_first, ctx) for level, s, is_first in level_plan(dict_config)]
    _validate_plan(plan)
    if ctx["no_compute"]:
        return []

    for p in plan: 
        epsilon = ctx["tol_res"]*np.sqrt(p["m"])#*2**(L-oneLvl)

        if p["is_first"]:
            y, Z, t_samples = get_sample_details(spde_model, wr_model, p["m"], p["d"], p["level"], ctx["L_first"], p["s"], sampling_fname)
            A, t_matrix = get_mtx(wr_model, p["J_s"], Z, p["d"], p["level"], p["s"], datamtx_fname)
            if log_freq:
                print("   Computing weights ...")
            w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), p["J_s"])
            if log_freq:
                print("   Weighted minimization ...")
            result = wr_model.method(A, y, w, p["s"], epsilon, ctx["nb_iter"], print_every=log_freq) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
            t_recovery = [result.tWC, result.tUser, result.tSys]
            lvl_by_lvl_result.append(CSPDEResult(p["J_s"], p["N"], p["s"], p["m"], p["d"], Z, y, 0, w, result, t_samples, t_matrix, t_recovery, p["t_J"]))
            print("\n\tRecovery time: {0} \t Building the Matrix: {1} \t Computing the samples: {2} \t Constructing polynomial set: {3} \n".format(t_recovery, t_matrix, t_samples, p["t_J"]))
        else:
            y, Z, t_samples = get_sample_details(spde_model, wr_model, p["m"], p["d"], p["level"], ctx["L_first"], p["s"], sampling_fname)
            A, t_matrix = get_mtx(wr_model, p["J_s"], Z, p["d"], p["level"], p["s"], datamtx_fname)
            if log_freq:
                print("   Computing weights ...")
            w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), p["J_s"])    
            # print(" Weights are {}".format(w) )
            if log_freq:
                print("   Weighted minimization ...")
            result = wr_model.method(A, y, w, p["s"], epsilon, ctx["nb_iter"], print_every=log_freq) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
            t_recovery = [result.tWC, result.tUser, result.tSys]
            # result = wr_model.method(A, y_new-y_old, w, p["s"], np.sqrt(p["m"]) *epsilon, ctx["nb_iter"]) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
            lvl_by_lvl_result.append(CSPDEResult(p["J_s"], p["N"], p["s"], p["m"], p["d"], Z, y, 0, w, result, t_samples, t_matrix, t_recovery, p["t_J"]))
            print("\n\tRecovery time: {0} \t Building the Matrix: {1} \t Computing the samples: {2} \t Constructing polynomial set: {3} \n".format(t_recovery, t_matrix, t_samples, p["t_J"]))

    # results = [_run_level(p, ctx) for p in plan]
    
    return lvl_by_lvl_result


def get_mtx(wr_model, J_s, Z, d, oneLvl, sl, datamtx_fname): 
    if (datamtx_fname is None):  
        # Create sampling matrix and weights
        print("   Creating sample operator ...")
        t_start = u.time_things()
        A = wr_model.operator.create(J_s, Z)
        t_matrix = u.time_things(t_start)
    else:
        # Create sampling matrix and weights
        mtx_file = datamtx_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '.npy'
        t_fname = datamtx_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '_time.npy'
        if os.path.isfile(t_fname):
            print("   Loading precomputed sample operator from {0} ...".format(mtx_file))
            A, t_matrix = wr_model.operator.load(mtx_file, t_fname)
        else: 
            print("   Creating sample operator ...")
            t_start = u.time_things()
            A = wr_model.operator.create(J_s, Z)
            A.save(mtx_file)
            t_matrix = u.time_things(t_start)
            np.save(t_fname, t_matrix)
            
    return A, t_matrix



def get_sample_details(spde_model, wr_model, m, d, oneLvl, J, sl, sampling_fname):
    if (sampling_fname is None):  
        Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
        print("\nComputing {0} SPDE sample approximations ...".format(m))
        # Get samples
        t_start = u.time_things()
        if oneLvl != J:
            y_old = spde_model.samples(Z)
            spde_model.refine_mesh()
        else:
            y_old = np.zeros(m)
        y_new = spde_model.samples(Z)
        t_samples = u.time_things(t_start)

    else:
        sampling_file = sampling_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '.npy'
        y_file = sampling_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '_y.npz'
        # sampling_file = sampling_fname + '_d' + str(d) + '_l' + str(oneLvl) + '.npy' ## Really HAVE to do this better one day!
        if os.path.isfile(sampling_file):
            Z = np.load(sampling_file)
        else: 
            Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
            np.save(sampling_file, Z)

        if os.path.isfile(y_file):
            print("\nLoading precomputed SPDE sample approximations from {0} ...".format(y_file))
            loaded = np.load(y_file)
            y_new = loaded['y_new']
            y_old = loaded['y_old']
            spde_model.refine_mesh()
            t_samples = [float(t) for t in loaded['t_samples']] if 't_samples' in loaded else [0.0, 0.0, 0.0]
        else:
            print("\nComputing {0} SPDE sample approximations ...".format(m))
            # Get samples
            t_start = u.time_things()
            if oneLvl != J:
                y_old = spde_model.samples(Z)
                spde_model.refine_mesh()
            else:
                y_old = np.zeros(m)
            y_new = spde_model.samples(Z)
            t_samples = u.time_things(t_start)
            np.savez(y_file, y_new=y_new, y_old=y_old, t_samples=t_samples)

    return y_new - y_old, Z, t_samples




def index_set(s, wr_model, ansatz_space):
    t = u.time_things()
    if ansatz_space == 0:
        J_s = weighted_sparse_indices(s, wr_model.operator.theta, wr_model.weights)
    else:
        dim = int(np.sum(np.isfinite(wr_model.weights)))
        J_s = total_deg_indices(dim, ansatz_space)
    t = u.time_things(t)
    return J_s, t

def total_deg_indices(v, max_degree = 2, threshold = np.inf):
    """All multi-indices nu in N^dim with |nu|_1 <= max_degree, without the (deg+1)^dim grid."""
    out = []
    nu = np.zeros(dim, dtype=int)

    def rec(k, budget):
        if k == dim:
            out.append(nu.copy())
            return
        for a in range(budget + 1):
            nu[k] = a
            rec(k + 1, budget - a)
        nu[k] = 0

    rec(0, max_degree)
    return out


def weighted_sparse_indices(s, theta, v):
    print("s = {0}, theta = {1}, v = {2}".format(s,theta,v))
    # Function for generating all admissible indices over given index set S
    def iterate(M, a, B, S):
        def iterate_(B, S, p):
            L = []

            if len(S):
                while True:
                    # Take the highest index in S
                    r     = S[-1]

                    # Substract weight from B
                    B    -= a[r]

                    # Increase nu_r by one
                    p[r] += 1

                    
                    if len(S) > 1:
                        # If there is more than one index left, recurse with remaining indices and 'remaining weight' B
                        e = iterate_(B, S[0:-1], p.copy())
                    else:
                        # If B - sum over a_j is larger than 0, add this multiindex
                        if B >= 0:
                            e = [p.copy()]
                        else:
                            e = []

                    # If the list of new multiindices is empty, there is nothing left to be done
                    # thanks to the monotonicity of a
                    if not e:
                        break

                    # Add found multiindices
                    L += e

            return L

        return iterate_(B, S, np.zeros(M, dtype='int'))


    # Set A and a as in Theorem 5.2
    A = np.log2(s/2.)
    a = 2 * np.log2(v)
    T = 2 * np.log2(theta)

    # Determine maximal M s.t. for j = 0 ... M-1 is a_j <= A - T
    # M is also the maximal support size
    # M = np.argmin(a <= A - T)
    M = np.argmin(a <= A )
    assert 0 != M, "Weight array too short. (Last element: {0}. Threshold: {1})".format(a[-1], A-T)

    # If A is non-negative the zero vector is always admissible
    assert A >= 0, "Negative A, i.e. sparsity less than 2."
    L = [np.zeros(M, dtype='int')]

    # Iterate through support sets of cardinality k = 1 ... M
    for k in range(1, M + 1):
        new_indices = []
    
        for S in itertools.combinations(range(M), k):
            new_indices += iterate(M, a, A - k*T, list(S))

        if [] == new_indices:
            break

        L += new_indices
    
    return L


def calculate_weights(theta, v, J_s):
    """omega_nu = theta^{|supp nu|} prod_j v_j^{nu_j}."""
    Ja = np.asarray(J_s)
    v = np.asarray(v, dtype=float)[:Ja.shape[1]]
    factors = np.where(Ja > 0, v[None, :] ** Ja, 1.0)
    return theta ** np.count_nonzero(Ja, axis=1) * np.prod(factors, axis=1)
