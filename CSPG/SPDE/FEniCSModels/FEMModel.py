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
__lastmodified__ = "2026/07/21"

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
        self.V = FunctionSpace(self.mesh, self.pde_cfg.get("elements", "Lagrange"), self.pde_cfg.get("degree", 1))

        # Define boundary conditions
        self.bc = DirichletBC(self.V, Constant(0.0), lambda x, on_boundary: on_boundary)

        # Define variational problem
        self.w = TrialFunction(self.V)
        self.v = TestFunction(self.V)

    def set_solver_parameters(self):
        """
        Parse the PDE configuration dictionary to set solver parameters.
        Parameters
        ----------
        pde_cfg : dict
            Dictionary containing PDE solver configuration parameters.

        TODO: Add more parameters to be parsed as needed.
        TODO: Propagate the use of pde config dict to other FEMModels.
        """
        self.linear_solver = self.pde_cfg.get("linear_solver", "gmres")
        self.preconditioner = self.pde_cfg.get("preconditioner", "amg")
        self.relative_tolerance = self.pde_cfg.get("relative_tolerance", 1e-6)
        self.absolute_tolerance = self.pde_cfg.get("absolute_tolerance", 1e-10)


    def refine_mesh(self, ratio=2): # Note, this can also be used to coarsen the mesh
        self.mesh_size = tuple(int(one_direction*ratio) for one_direction in self.mesh_size)
        self.init_simple_mesh()

    def set_mesh_size(self, mesh_size):
        self.mesh_size = tuple(one_size for one_size in mesh_size)
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
        # print assemble(self.M)
        # return np.sum(u.vector().get_local())
        return assemble(self.M)
        # return self.solver.evaluate_goal(Form(self.M), u)

    def __getstate__(self):
        odict = self.__dict__.copy()

        # Can't pickle dolfin's mesh, solver and M...
        for key in ['mesh', 'solver', 'M', 'V', 'bc', 'w', 'v']:
            if key in odict:
                del odict[key]
        return odict
