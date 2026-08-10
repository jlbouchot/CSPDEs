from dolfin import *
from .FEMModel import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2017/10/25"

class DiffusionFEMModelML(FEMModel):
    def __init__(self, a, f, M_gen, mesh_size, pde_cfg):
        self.a         = a
        self.f         = f
        self.M_gen     = M_gen
        self.pde_cfg    = pde_cfg

        self.mesh_size = mesh_size
        self.init_simple_mesh()
        # self.parse_pde_cfg(pde_cfg)

    def solve(self, z):
        # TODO: Make this a bit more elegant, likely letting the log level being specified on the fly, and moving the various values to another file.
        # Make FEniCS output only the most important messages
        CRITICAL  = 50 #, // errors that may lead to data corruption and suchlike
        ERROR     = 40 #, // things that go boom
        WARNING   = 30 #, // things that may go boom later
        INFO      = 20 #, // information of general interest
        PROGRESS  = 16 #, // what's happening (broadly)
        TRACE     = 13 #, // what's happening (in detail)
        DBG       = 10#   // sundry
        set_log_level(ERROR)

        # Create mesh if there is none
        if not hasattr(self, 'mesh'):
            self.init_simple_mesh()

        params = self.split_params([self.a, self.f], z)

        x = SpatialCoordinate(self.mesh)
        A = self.a(x, Constant(params[0])) * inner(nabla_grad(self.w), nabla_grad(self.v)) * dx
        L = self.f(x, Constant(params[1])) * self.v * dx

        # Create goal-functional for error estimation
        u      = Function(self.V)
        self.M = self.M_gen(self, u, dx)

        # Create solver
        problem     = LinearVariationalProblem(A, L, u, self.bc)
        self.solver = LinearVariationalSolver(problem)
        # self.solver = LinearVariationalSolver(problem, solver_parameters={
        #   "linear_solver": "gmres", 
        #   "preconditioner": "amg"})
        self.set_solver_parameters()
        # self.solver.parameters["linear_solver"] = self.linear_solver
        # self.solver.parameters["preconditioner"] = self.preconditioner
        # self.solver.parameters.add("relative_tolerance", self.relative_tolerance)
        # self.solver.parameters.add("absolute_tolerance", self.absolute_tolerance)

        # Compute solution
        self.solver.solve()

        return u
