from dolfin import *
import numpy as np
from ..     import SPDEModel

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

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
                assert False, "Only one to three dimensional problems supported"
        else:
            self.mesh = UnitIntervalMesh(self.mesh_size)

    def refine_mesh(self, ratio=2): # Note, this can also be used to coarsen the mesh
        self.mesh_size = tuple(int(one_direction*ratio) for one_direction in self.mesh_size)
        self.init_simple_mesh()

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
        return assemble(self.M)
        # return self.solver.evaluate_goal(Form(self.M), u)

    def __getstate__(self):
        odict = self.__dict__.copy()

        # Can't pickle mesh, solver and M
        del odict['mesh']
        del odict['solver']
        del odict['M']

        return odict
