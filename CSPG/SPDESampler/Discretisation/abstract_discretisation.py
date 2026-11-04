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

from abc import ABC, abstractmethod

class AbstractDiscretisation(ABC):
    '''Abstract class for a discretisation of a PDE. It contains the mesh, the finite element space, the Dirichlet boundary conditions, and the free dofs. It also provides methods to assemble matrices and vectors, and to convert between free dofs and full functions.'''

    def __init__(self, mesh_size, degree=1):
        self.mesh_size = None
        self.mesh = None
        self.V = None
        self.bc = None
        self.x = None
        self.free = None
        self._gather = None
        self.is_unit = False

    def matrix_values(self, form):
        '''Assemble a bilinear form, return the CSR values on the free dofs (pattern fixed at the first call).'''
        pass

    def vector_values(self, form):
        '''Assemble a linear form, return the values on the free dofs.'''
        pass

    def function(self, u_free):
        '''Convert a vector of free dofs to a full function in the finite element space.'''
        pass
