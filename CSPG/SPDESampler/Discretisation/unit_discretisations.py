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
import scipy.sparse as sp
from dolfin import *

from .abstract_discretisation import AbstractDiscretisation

class UnitDiscretisation(AbstractDiscretisation):
    '''Mesh + P_k space + elimination of the (homogeneous) Dirichlet dofs. The reduced matrices stay SPD.'''

    def __init__(self, mesh_size, degree=1):
        mesh_cls = {1: UnitIntervalMesh, 2: UnitSquareMesh, 3: UnitCubeMesh}[len(mesh_size)]
        self.mesh_size = tuple(mesh_size)
        self.mesh = mesh_cls(*mesh_size)
        self.V = FunctionSpace(self.mesh, "Lagrange", degree)
        self.bc = DirichletBC(self.V, Constant(0.0), lambda x, on_boundary: on_boundary)
        self.x = SpatialCoordinate(self.mesh)
        bnd = np.fromiter(self.bc.get_boundary_values().keys(), dtype=np.int64)
        self.free = np.setdiff1d(np.arange(self.V.dim()), bnd)
        self._gather = None
        self.is_unit = True

    def matrix_values(self, form):
        '''Assemble a bilinear form, return the CSR values on the free dofs (pattern fixed at the first call).'''
        ip, ix, data = as_backend_type(assemble(form)).mat().getValuesCSR()
        if self._gather is None:
            n = ip.size - 1
            pos = sp.csr_matrix((np.arange(1, ix.size + 1, dtype=float), ix, ip), shape=(n, n))[self.free][:, self.free].tocsr()
            pos.sort_indices()
            self._gather, self.indptr, self.indices = pos.data.astype(np.int64) - 1, pos.indptr, pos.indices
        return data[self._gather]

    def vector_values(self, form):
        return assemble(form).get_local()[self.free]

    def function(self, u_free):
        u = Function(self.V)
        vec = np.zeros(self.V.dim())
        vec[self.free] = u_free
        u.vector().set_local(vec)
        u.vector().apply("insert")
        return u
