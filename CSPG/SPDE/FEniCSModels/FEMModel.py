from dolfin import *
import numpy as np
from ..     import SPDEModel

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/01/22"

class FEMModel(SPDEModel):
    def init_simple_mesh(self):
        # Create mesh and define function space
        if type(self.mesh_size) is tuple:
            if len(self.mesh_size) == 1:
                self.mesh = UnitIntervalMesh(*self.mesh_size)
            elif len(self.mesh_size) == 2:
                self.mesh = UnitSquareMesh(*self.mesh_size)
            elif len(self.mesh_size) == 3:
                self.mesh = UnitCubeMesh(*self.mesh_size)
            else:
                assert False, f"Init simple mesh passed a {len(self.mesh_size)} dimensional problem. Only 1 to 3 dimensions supported"
        else:
            self.mesh = UnitIntervalMesh(self.mesh_size)

        self.generate_functions_spaces()

    def generate_functions_spaces(self): 
        # Create approximation space
        self.V = FunctionSpace(self.mesh, 'Lagrange', 2)

        # Define boundary conditions
        self.bc = DirichletBC(self.V, Constant(0.0), lambda x, on_boundary: on_boundary)

        # Define variational problem
        self.w = TrialFunction(self.V)
        self.v = TestFunction(self.V)


    def refine_mesh(self, ratio=2): # Note, this can also be used to coarsen the mesh
        self.mesh_size = tuple(int(one_direction*ratio) for one_direction in self.mesh_size)
        self.init_simple_mesh()

    def set_mesh_size(self, mesh_size):
        self.mesh_size = tuple(one_size for one_size in mesh_size)
        self.init_simple_mesh

    # @staticmethod
    def split_params(self, coeff, z):
        cur_split = 0
        params    = []

        for c in coeff:
            z_c = np.array(0)

            if c.num_params > 0:
                z_c = z[cur_split:(cur_split+c.num_params)]
                cur_split += c.num_params

            params.append(z_c)
        return params

    def sample(self, z):
        u = self.solve(z)
        # Return functional value for solution
        # print assemble(self.M)
        # return np.sum(u.vector().get_local())
        return assemble(self.M)
        # return self.solver.evaluate_goal(Form(self.M), u)

    def __getstate__(self):
        odict = self.__dict__.copy()

        # Can't pickle dolfin's mesh, solver and M...
        del odict['mesh']
        del odict['solver']
        del odict['M']
        del odict['V']
        del odict['bc']
        del odict['w']
        del odict['v']

        return odict
