"""Shared functionality for PDE samplers.

The base class owns the discretisation, parameter-dependent coefficient, and
quantity-of-interest evaluation. Concrete samplers implement ``solve_free``
to solve their discrete PDE and return the free degrees of freedom.
"""
__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/29"
__lastmodified__ = "2026/11/04"

import numpy as np
from dolfin import Constant, TestFunction, assemble
from ..Solvers.linear_solvers import make_solver
from progressbar import Bar, ETA, Percentage, ProgressBar

class BaseSampler:
    """Base class for samplers that solve a parameterized PDE."""

    def __init__(self, disc, physics, qoi, solver):
        """Initialize with a discretisation, PDE physics, QoI, and solver."""
        self.disc, self.physics, self.qoi = disc, physics, qoi
        self.solver = make_solver(solver) if isinstance(solver, str) else solver
        self.zc = Constant(np.zeros(max(physics.num_params, 1)))
        self.c = physics.coefficient(disc.x, self.zc)
        self.g = disc.vector_values(qoi.form(disc, TestFunction(disc.V))) if qoi.linear else None
        self.stats = dict(newton_steps=0)

    def set_z(self, z):
        """Set the active parameter vector, padding or truncating to its dimension."""
        zz = np.zeros(max(self.physics.num_params, 1))
        k = min(zz.size, len(z))
        zz[:k] = z[:k]                      # extra parameters are inert (as the legacy split_params)
        self.zc.assign(Constant(zz))
        return zz

    def value(self, u_free):
        """Evaluate the quantity of interest for a solution in free-dof form."""
        if self.g is not None:
            return float(self.g @ u_free)
        return assemble(self.qoi.form(self.disc, self.disc.function(u_free)))

    def solve_free(self, z):
        """Solve the PDE for ``z`` and return solution values on free dofs."""
        raise NotImplementedError

    def sample(self, z):
        """Solve the PDE and return its quantity-of-interest value."""
        return self.value(self.solve_free(z))

    def samples(self, Z):

        m, _ = Z.shape

        # Show a progressbar
        widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
        pbar    = ProgressBar(widgets=widgets)

        # Compute all functional evalutations with given precision
        y = np.zeros(m)
        for k in pbar(range(m)):
            y[k] = self.sample(Z[k])

        return y
