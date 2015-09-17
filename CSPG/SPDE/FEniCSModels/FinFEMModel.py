from dolfin import *
from FileFEMModel import *

import os


class FinFEMModel(FileFEMModel):
    def __init__(self, a, Bi, M_gen):
        self.a  = a
        self.Bi = Bi

        self.M_gen = M_gen

        self.init_mesh("fin")

    def solve(self, z, epsilon):
        # Make FEniCS output only the most important messages
        set_log_level(WARNING)

        # Create mesh if necessary and define function space
        # if not hasattr(self, 'mesh'):
        #     self.init_mesh()

        # TODO: Reusing adapted mesh does not work here.
        #       (Crash in next call to solve if mesh has been adapted.)
        self.init_mesh("fin")

        V = FunctionSpace(self.mesh, "Lagrange", 2)

        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)

        dx = Measure("dx")[self.subdomains]
        ds = Measure("ds")[self.boundaries]

        x      = SpatialCoordinate(self.mesh)
        params = self.split_params(self.a + [self.Bi], z)

        a = self.a[0](x, Constant(params[0])) * inner(grad(u), grad(v)) * dx(7)   + \
            self.a[1](x, Constant(params[1])) * inner(grad(u), grad(v)) * dx(8)   + \
            self.a[2](x, Constant(params[2])) * inner(grad(u), grad(v)) * dx(9)   + \
            self.a[3](x, Constant(params[3])) * inner(grad(u), grad(v)) * dx(10)  + \
            self.a[4](x, Constant(params[4])) * inner(grad(u), grad(v)) * dx(11)  + \
            self.Bi  (x, Constant(params[5])) * u * v * (ds(2) + ds(3) + ds(4) + ds(5) + ds(6))
        L = v * ds(1)
        F = a - L

        # Separate left and right hand sides of equation
        a, L = lhs(F), rhs(F)

        # Compute solution
        u = Function(V)

        self.M      = self.M_gen(self, u, ds(1))
        problem     = LinearVariationalProblem(a, L, u)
        self.solver = AdaptiveLinearVariationalSolver(problem, self.M)
        self.solver.solve(epsilon)

        return u
