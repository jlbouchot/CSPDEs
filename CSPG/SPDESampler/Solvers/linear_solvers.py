from abc import ABC, abstractmethod

import numpy as np
import utilities.utils as utils

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/29"
__lastmodified__ = "2026/11/07"

REGISTRY = {}
# dolfin "linear_solver" names used in the .ini files -> registry names
ALIASES = {"default": "umfpack", "lu": "umfpack", "petsc": "petsc-lu", "superlu_dist": "superlu"}


def register(cls):
    REGISTRY[cls.name] = cls
    return cls


def make_solver(name, **options):
    name = ALIASES.get(str(name).lower(), str(name).lower())
    if name not in REGISTRY:
        raise ValueError("Unknown linear solver '{0}'. Available: {1}".format(name, sorted(REGISTRY)))
    return REGISTRY[name](**options)


class LinearSolver(ABC):
    name = None
    exact = True              # False for iterative solvers: results depend on a tolerance
    requires_spd = False      # True: give it the SPD matrix (Dirichlet dofs eliminated, not "bc.apply"-ed)
    reuses_symbolic = True    # False: every factorize() redoes the ordering

    def __init__(self, **options):
        self.options = options
        self.indptr = self.indices = None
        self.times = {}

    def set_pattern(self, indptr, indices):
        self.reset()
        self.indptr = np.asarray(indptr, dtype=np.int32)
        self.indices = np.asarray(indices, dtype=np.int32)
        self.n = self.indptr.size - 1
        t = utils.time_things()
        self._analyze()
        self.times["analyze"] = utils.time_things(t)

    def factorize(self, data):
        assert self.indptr is not None, "set_pattern() first"
        data = np.asarray(data, dtype=np.float64)
        assert data.size == self.indices.size, "values do not match the pattern"
        t = utils.time_things()
        self._factorize(data)
        self.times["factorize"] = utils.time_things(t)

    def solve(self, rhs):
        t = utils.time_things()
        u = self._solve(np.asarray(rhs, dtype=np.float64))
        self.times["solve"] = utils.time_things(t)
        # return (u, self.times) # TODO: Eventually return detailed timing information. For now, just return the solution vector.
        return u
        
    def reset(self):
        """Forget the pattern (mesh change)."""
        self.indptr = self.indices = None

    # Needed for any concrete solver implementation. 
    def _analyze(self):
        pass

    @abstractmethod
    def _factorize(self, data):
        ...

    @abstractmethod
    def _solve(self, rhs):
        ...

    # Some helper functions

    def csr(self, data):
        import scipy.sparse as sp
        return sp.csr_matrix((data, self.indices, self.indptr), shape=(self.n, self.n))

    def __repr__(self):
        return "{0}({1})".format(type(self).__name__, self.options or "")


# ====================================================================== PETSc back-ends (petsc4py)
class _PETScFactor(LinearSolver):
    """KSPPREONLY + PCLU/PCCHOLESKY: PETSc keeps ordering + symbolic factorisation while the Mat pattern is unchanged."""
    package = None
    factor = "lu"
    ordering = None

    def _analyze(self):
        from petsc4py import PETSc
        self._PETSc = PETSc
        self._mat = None

    def _factorize(self, data):
        PETSc = self._PETSc
        if self._mat is None:
            self._mat = PETSc.Mat().createAIJ((self.n, self.n), csr=(self.indptr, self.indices, data), comm=PETSc.COMM_SELF)
            self._mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)
            if self.requires_spd:
                self._mat.setOption(PETSc.Mat.Option.SPD, True)
            self._mat.assemble()
            ksp = PETSc.KSP().create(PETSc.COMM_SELF)
            ksp.setType("preonly")
            pc = ksp.getPC()
            pc.setType(self.factor)
            pc.setFactorSolverType(self.package)
            if self.ordering:
                pc.setFactorOrdering(self.ordering)
            ksp.setOperators(self._mat)
            self._ksp = ksp
            self._b, self._x = self._mat.createVecLeft(), self._mat.createVecRight()
        else:
            self._mat.setValuesCSR(self.indptr, self.indices, data)
            self._mat.assemble()
        self._ksp.setUp()                      # numeric factorisation happens here (symbolic reused)

    def _solve(self, rhs):
        self._b.setArray(rhs)
        self._ksp.solve(self._b, self._x)
        if self._ksp.getConvergedReason() < 0:
            raise RuntimeError("{0}: solve failed (reason {1})".format(self.name, self._ksp.getConvergedReason()))
        return self._x.getArray().copy()

    def reset(self):
        super().reset()
        if getattr(self, "_mat", None) is not None:
            self._ksp.destroy()
            self._mat.destroy()
        self._mat = None


@register
class UmfpackLU(_PETScFactor):
    name, package = "umfpack", "umfpack"


@register
class MumpsLU(_PETScFactor):
    name, package = "mumps", "mumps"          # do not use in fork()ed workers (MPI): use spawn


@register
class MumpsCholesky(_PETScFactor):
    name, package, factor, requires_spd = "mumps-cholesky", "mumps", "cholesky", True


@register
class PETScLU(_PETScFactor):
    name, package, ordering = "petsc-lu", "petsc", "nd"


@register
class SuperLU(_PETScFactor):
    name, package = "superlu", "superlu"


@register
class CholmodPETSc(_PETScFactor):
    """CHOLMOD through PETSc (SuiteSparse build): the fastest serial solver measured on these systems."""
    name, package, factor, requires_spd = "cholmod", "cholmod", "cholesky", True


@register
class CGAMG(LinearSolver):
    """Iterative: CG + hypre BoomerAMG (falls back to GAMG). Not exact: only for the ground truth / very fine meshes."""
    name, exact, requires_spd = "cg-amg", False, True

    def _analyze(self):
        from petsc4py import PETSc
        self._PETSc, self._mat = PETSc, None
        ## How often the AMG hierarchy is rebuilt, counted in factorize() calls:
        ##   1 (default) every sample -- what PETSc does if left alone
        ##   k > 1       every k-th sample, reusing it in between
        ##   0           once, then never again
        ## The hierarchy setup dominates an AMG solve, and here the sparsity
        ## pattern never changes between samples, only the coefficient values,
        ## so reusing it is the main speed lever. The price is a preconditioner
        ## that is only approximate for the new matrix, i.e. more CG iterations
        ## and, if pushed too far, a solve that no longer reaches rtol. Do not
        ## trust it without measuring: compute_ground_truth.py --calibrate
        ## reports the iteration count AND the QoI error against a direct solve.
        self._rebuild_every = int(self.options.get("pc_rebuild_every", 1))
        self._nfac = 0

    def _factorize(self, data):
        PETSc = self._PETSc
        if self._mat is None:
            self._mat = PETSc.Mat().createAIJ((self.n, self.n), csr=(self.indptr, self.indices, data), comm=PETSc.COMM_SELF)
            self._ksp = PETSc.KSP().create(PETSc.COMM_SELF)
            self._ksp.setType("cg")
            pc = self._ksp.getPC()
            try:
                pc.setType("hypre")
            except PETSc.Error:
                pc.setType("gamg")
            self._ksp.setTolerances(rtol=self.options.get("rtol", 1e-10), atol=self.options.get("atol", 1e-14))
            self._ksp.setOperators(self._mat)
            self._b, self._x = self._mat.createVecLeft(), self._mat.createVecRight()
        else:
            self._mat.setValuesCSR(self.indptr, self.indices, data)
        self._mat.assemble()
        ## Left completely alone when _rebuild_every == 1, so the default path
        ## through this solver is byte for byte what it was.
        if self._rebuild_every != 1:
            fresh = self._nfac == 0 or (self._rebuild_every > 0
                                        and self._nfac % self._rebuild_every == 0)
            self._ksp.getPC().setReusePreconditioner(not fresh)
        self._nfac += 1

    def _solve(self, rhs):
        self._b.setArray(rhs)
        self._ksp.solve(self._b, self._x)
        return self._x.getArray().copy()


# ====================================================================== pure SciPy (no PETSc)
@register
class ScipySuperLU(LinearSolver):
    """scipy.sparse.linalg.splu. No symbolic reuse, but no PETSc either: safe in any worker process."""
    name, reuses_symbolic = "scipy", False

    def _factorize(self, data):
        import scipy.sparse.linalg as spla
        self._lu = spla.splu(self.csr(data).tocsc(), permc_spec="MMD_AT_PLUS_A")

    def _solve(self, rhs):
        return self._lu.solve(rhs)


# ====================================================================== hand-made example
@register
class CholmodSkSparse(LinearSolver):
    """Template for a hand-made solver: CHOLMOD called directly through scikit-sparse.

    analyze() once per pattern (AMD/METIS ordering + supernodal symbolic factorisation),
    cholesky_inplace() per sample (numeric only). `pip install scikit-sparse` (needs SuiteSparse headers).
    """
    name, requires_spd = "cholmod-sksparse", True

    def _analyze(self):
        self._factor = None

    def _factorize(self, data):
        A = self.csr(data).tocsc()
        if self._factor is None:
            from sksparse.cholmod import analyze
            self._factor = analyze(A, mode=self.options.get("mode", "supernodal"), ordering_method=self.options.get("ordering", "default"))
        self._factor.cholesky_inplace(A)

    def _solve(self, rhs):
        return self._factor(rhs)
