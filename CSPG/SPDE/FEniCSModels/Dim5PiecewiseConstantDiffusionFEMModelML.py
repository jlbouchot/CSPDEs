from dolfin import *
from .FEMModel import *
import numpy as np

# Still have to find the way to 
# 1) Automatize the process to an unknown d number of sub domains and 
# 2) uncertain splitting points
class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 0 and x <= 1.0/5.0 else False

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 1.0/5.0 and x <= 2.0/5.0 else False

class Omega2(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 2.0/5.0 and x <= 3.0/5.0 else False

class Omega3(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 3.0/5.0 and x <= 4.0/5.0 else False

class Omega4(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 4.0/5.0 and x <= 5.0/5.0 else False


tol = 1E-14   # tolerance for coordinate comparisons
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x) < tol


class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x - 1) < tol


class Dim5PiecewiseConstantDiffusionFEMModelML(FEMModel):
    
    def __init__(self, abar, f, var, mesh_size, M_gen):
        
        self.abar      = abar
        self.f         = f
        self.var       = var
        self.M_gen     = M_gen
        self.num_subd  = 6

        self.mesh_size = mesh_size
        self.init_simple_mesh()
        # I'm assuming, this is where we have to define the different regions... kinda annoying.

    def solve(self, z):
        # Make FEniCS output only the most important messages
        set_log_level(WARNING)

        # Create mesh if there is none
        if not hasattr(self, 'mesh'):
            self.init_simple_mesh()

        # Create the d subdomains
        subdomains = MeshFunction('size_t', self.mesh, 1) # The last argument corresponds to the dimension of the cells: here, intervals, dim 1. 
        subdomains.set_all(0) ### IMPORTANT! This is not a GOOD way to deal with it. But I have no clue how to solve the problem. It appears that the 'meshes' at the boundaries between two partitions are marked with random numbers
        subdomain0 = Omega0()
        subdomain0.mark(subdomains, 0)
        subdomain1 = Omega1()
        subdomain1.mark(subdomains, 1)
        subdomain2 = Omega2()
        subdomain2.mark(subdomains, 2)
        subdomain3 = Omega3()
        subdomain3.mark(subdomains, 3)
        subdomain4 = Omega4()
        subdomain4.mark(subdomains, 4)

        V0 = FunctionSpace(self.mesh, 'DG', 0) # Function space of constant functions
        # This has to be used for the uncertain diffusion coefficients
        k  = Function(V0) # That particular -- parametric -- diffusion coefficient
        k_values = z*self.var+self.abar
        # Affect the appropriate local diffusion value:
        help = np.asarray(subdomains.array(), dtype=np.int32)
        k.vector()[:] = np.choose(help, k_values)		

        # Now we can keep going with the usual FEniCS process.

        # Create approximation space
        V = FunctionSpace(self.mesh, 'Lagrange', 2)
        Gamma_0 = DirichletBC(V, Constant(0), BottomBoundary())
        Gamma_1 = DirichletBC(V, Constant(1), TopBoundary())

        # Define boundary conditions
        bc = [Gamma_0, Gamma_1] # DirichletBC(V, Constant(0.0), lambda x, on_boundary: on_boundary)

        # Define variational problem
        w = TrialFunction(V)
        v = TestFunction(V)


        x = SpatialCoordinate(self.mesh)
        A = k * inner(nabla_grad(w), nabla_grad(v)) * dx
        L = self.f(x, Constant(1)) * v * dx

        # Create goal-functional for error estimation
        u      = Function(V)
        self.M = self.M_gen(self, u, dx)

        # Create solver
        problem     = LinearVariationalProblem(A, L, u, bc)
        self.solver = LinearVariationalSolver(problem)
        # y[k] = assemble(myAverage(mesh, u, dx))

        # Compute solution
        self.solver.solve()

        return u
