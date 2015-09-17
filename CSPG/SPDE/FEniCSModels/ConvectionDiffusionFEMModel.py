from dolfin import *
from FileFEMModel import *

import os


class ConvectionDiffusionFEMModel(FileFEMModel):
    def __init__(self, a, b, M_gen):
        self.a = a
        self.b = b

        self.M_gen = M_gen

        self.init_mesh("fun_square")

    def solve(self, z, epsilon):
        # Make FEniCS output only the most important messages
        set_log_level(WARNING)

        # Create mesh if necessary and define function space
        # if not hasattr(self, 'mesh'):
        #     self.init_mesh()

        # TODO: Reusing adapted mesh does not work here.
        #       (Crash in next call to solve if mesh has been adapted.)
        # self.init_mesh("fun_square")
        self.init_mesh("cv_square")

        V = FunctionSpace(self.mesh, "Lagrange", 3)

        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)

        dx = Measure("dx")[self.subdomains]
        ds = Measure("ds")[self.boundaries]

        x      = SpatialCoordinate(self.mesh)
        params = self.split_params([self.a, self.b], z)

        bcs = DirichletBC(V, Constant(0.0), self.boundaries, 1)

        f = Expression("20 * exp(-(pow(x[0] - 1.5, 2) + pow(x[1] - 1.5, 2)) / 0.5)")
        g = Constant(0.0)
        # TODO: "b" is not supposed to be "Constant"
        a = self.a(x, Constant(params[0])) * inner(grad(u), grad(v)) * (dx(3)+dx(4)) \
            + inner(Constant(self.b(x, params[1])), grad(u)) * v     * (dx(3)+dx(4))
        L = (f * v * (dx(3)+dx(4))) + g * v * ds(2)
        F = a - L

        # Separate left and right hand sides of equation
        a, L = lhs(F), rhs(F)

        # Compute solution
        u = Function(V)

        self.M      = self.M_gen(self, u, dx(4))
        problem     = LinearVariationalProblem(a, L, u, bcs)
        self.solver = AdaptiveLinearVariationalSolver(problem, self.M)
        self.solver.solve(epsilon)

        return u
