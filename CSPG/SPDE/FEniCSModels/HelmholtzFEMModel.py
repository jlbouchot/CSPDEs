from dolfin import *
from .FileFEMModel import *


class HelmholtzFEMModel(FileFEMModel):
    def __init__(self, rho, k, M_gen):
        self.rho = rho
        self.k   = k

        self.M_gen = M_gen

        self.init_mesh("square")

    def solve(self, z, epsilon):
        # Make FEniCS output only the most important messages
        set_log_level(WARNING)

        # Create mesh if necessary and define function space
        # if not hasattr(self, 'mesh'):
        #     self.init_mesh()

        # TODO: Reusing adapted mesh does not work here.
        #       (Crash in next call to solve if mesh has been adapted.)
        self.init_mesh("square")

        # Create function spaces
        V = FunctionSpace(self.mesh, 'CG', 4)
        W = V*V

        # Create source and boundary conditions
        f_r, f_i   = Expression("exp(-(pow(x[0] - 2.5, 2) + pow(x[1] - 2.5, 2)) / 0.01)"), Constant(0.0)
        u0_r, u0_i = Constant(0.0), Constant(0.0)

        bcs = [DirichletBC(W.sub(0), u0_r, self.boundaries, 1),
               DirichletBC(W.sub(1), u0_i, self.boundaries, 1)]

        # Define variational problem
        dx = Measure("dx")[self.subdomains]
        ds = Measure("ds")[self.boundaries]

        (u_r, u_i) = TrialFunctions(W)
        (v_r, v_i) = TestFunctions(W)

        x      = SpatialCoordinate(self.mesh)
        params = self.split_params([self.rho] + self.k, z)

        a_r = inner(grad(u_r), grad(v_r)) * 1./self.rho   (x, Constant(params[0]))     * (dx(3)+dx(4)) \
            - inner(grad(u_i), grad(v_i)) * 1./self.rho   (x, Constant(params[0]))     * (dx(3)+dx(4)) \
            - (u_r * v_r - u_i * v_i)     *    self.k  [0](x, Constant(params[1]))**2  * dx(3) \
            - (u_r * v_r - u_i * v_i)     *    self.k  [1](x, Constant(params[2]))**2  * dx(4) \
            + inner(grad(u_r), grad(v_r)) * 1./self.rho   (x, Constant(params[0]))     * ds(1) \
            - inner(grad(u_i), grad(v_i)) * 1./self.rho   (x, Constant(params[0]))     * ds(1) \
            - (u_r * v_i + u_i * v_r)     *    self.k  [0](x, Constant(params[1]))**2  * ds(1)
            

        a_i = inner(grad(u_r), grad(v_i)) * 1./self.rho   (x, Constant(params[0]))    * (dx(3)+dx(4)) \
            + inner(grad(u_i), grad(v_r)) * 1./self.rho   (x, Constant(params[0]))    * (dx(3)+dx(4)) \
            - (u_r * v_i + u_i * v_r)     *    self.k  [0](x, Constant(params[1]))**2 * dx(3) \
            - (u_r * v_i + u_i * v_r)     *    self.k  [1](x, Constant(params[2]))**2 * dx(4) \
            + inner(grad(u_r), grad(v_i)) * 1./self.rho   (x, Constant(params[0]))    * ds(1) \
            + inner(grad(u_i), grad(v_r)) * 1./self.rho   (x, Constant(params[0]))    * ds(1) \
            + (u_r * v_r - u_i * v_i)     *    self.k  [0](x, Constant(params[1]))**2 * ds(1)


        L_r = f_r * v_r * (dx(3)+dx(4)) - f_i * v_i * (dx(3)+dx(4))
        L_i = f_r * v_i * (dx(3)+dx(4)) + f_i * v_r * (dx(3)+dx(4))

        a = a_r + a_i
        L = L_r + L_i

        F = a - L
        a, L = lhs(F), rhs(F)

        # Compute solution
        u = Function(W)

        self.M      = self.M_gen(self, split(u)[0], dx(4))
        problem     = LinearVariationalProblem(a, L, u, bcs)
        self.solver = AdaptiveLinearVariationalSolver(problem, self.M)
        self.solver.solve(epsilon)

        return u
