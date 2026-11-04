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
from dolfin import Constant

from .base_sampler import BaseSampler

class AssembledLinearSampler(BaseSampler):
    '''Any linear problem: assemble A(z) per sample (pattern fixed => the solver keeps its symbolic factorisation).'''

    def __init__(self, *args):
        super().__init__(*args)
        self.b = self.disc.vector_values(self.physics.linear(self.disc))
        self.A_form = self.physics.bilinear(self.disc, self.c)
        self.disc.matrix_values(self.A_form)
        self.solver.set_pattern(self.disc.indptr, self.disc.indices)

    def solve_free(self, z):
        self.set_z(z)
        self.solver.factorize(self.disc.matrix_values(self.A_form))
        return self.solver.solve(self.b)


class AffineLinearSampler(AssembledLinearSampler):
    '''A(z) = A_0 + sum_k z_k A_k: p+1 assemblies per mesh, then one dense mat-vec on the CSR values per sample.'''

    def __init__(self, *args):
        super().__init__(*args)
        p = self.physics.num_params
        e = np.zeros(max(p, 1))
        self.zc.assign(Constant(e))
        self.D0 = self.disc.matrix_values(self.A_form)
        self.D = np.empty((p, self.D0.size))
        for k in range(p):
            e[:] = 0.0
            e[k] = 1.0
            self.zc.assign(Constant(e))
            self.D[k] = self.disc.matrix_values(self.A_form) - self.D0
        z = np.random.default_rng(0).uniform(-1, 1, p)           # check the affine claim once per mesh
        self.zc.assign(Constant(z))
        exact = self.disc.matrix_values(self.A_form)
        if np.linalg.norm(exact - self.D0 - z @ self.D) > 1e-10 * np.linalg.norm(exact):
            raise ValueError("coefficient is not affine in z: use AssembledLinearSampler")

    def solve_free(self, z):
        zz = np.zeros(self.physics.num_params)
        k = min(zz.size, len(z))
        zz[:k] = z[:k]
        self.solver.factorize(self.D0 + zz @ self.D)
        return self.solver.solve(self.b)
