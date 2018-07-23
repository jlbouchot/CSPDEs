from dolfin import *
from .FEMModel import *
import numpy as np

# Still have to find the way to 
# 1) Automatize the process to an unknown d number of sub domains and 
# 2) uncertain splitting points
class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 0 and x <= 1.0/13.0 else False

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 1.0/13.0 and x <= 2.0/13.0 else False

class Omega2(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 2.0/13.0 and x <= 3.0/13.0 else False

class Omega3(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 3.0/13.0 and x <= 4.0/13.0 else False

class Omega4(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 4.0/13.0 and x <= 5.0/13.0 else False

class Omega5(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 5.0/13.0 and x <= 6.0/13.0 else False

class Omega6(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 6.0/13.0 and x <= 7.0/13.0 else False

class Omega7(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 7.0/13.0 and x <= 8.0/13.0 else False

class Omega8(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 8.0/13.0 and x <= 9.0/13.0 else False

class Omega9(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 9.0/13.0 and x <= 10.0/13.0 else False

class Omega10(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 10.0/13.0 and x <= 11.0/13.0 else False

class Omega11(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 11.0/13.0 and x <= 12.0/13.0 else False

class Omega12(SubDomain):
    def inside(self, x, on_boundary):
        return True if x >= 12.0/13.0 and x <= 1.0 else False

tol = 1E-14   # tolerance for coordinate comparisons
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x) < tol


class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x - 1) < tol


class PiecewiseConstantDiffusionFEMModelML(FEMModel):
    
    def __init__(self, abar, f, var, mesh_size, M_gen):
        
        self.abar      = abar
        self.f         = f
        self.var       = var
        self.M_gen     = M_gen
        self.num_subd  = 13

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
        subdomain5 = Omega5()
        subdomain5.mark(subdomains, 5)
        subdomain6 = Omega6()
        subdomain6.mark(subdomains, 6)
        subdomain7 = Omega7()
        subdomain7.mark(subdomains, 7)
        subdomain8 = Omega8()
        subdomain8.mark(subdomains, 8)
        subdomain9 = Omega9()
        subdomain9.mark(subdomains, 9)
        subdomain10 = Omega10()
        subdomain10.mark(subdomains, 10)
        subdomain11 = Omega11()
        subdomain11.mark(subdomains, 11)
        subdomain12 = Omega12()
        subdomain12.mark(subdomains, 12)

        V0 = FunctionSpace(self.mesh, 'DG', 0) # Function space of constant functions
        # This has to be used for the uncertain diffusion coefficients
        k  = Function(V0) # That particular -- parametric -- diffusion coefficient
        ### Have to improve in the following part for the case of "more" uncertainty
        # # Recover the diffusion values from the parameter vector z passed as an argument:
        # params = self.split_params(self.a + [self.f], z)
        # # Hopefully params[0:8] should contain the coefs we are looking for 
        # k_values = params[0:8]/2+5 # The 5 should be changed to the abar given as a parameter
        k_values = z*self.var+self.abar
        # Affect the appropriate local diffusion value:
        help = np.asarray(subdomains.array(), dtype=np.int32)
        k.vector()[:] = np.choose(help, k_values)		

        # Now we can keep going with the usual FEniCS process.

        # Create approximation space
        V = FunctionSpace(self.mesh, 'Lagrange', 2)
        # V = FunctionSpace(self.mesh, 'Lagrange', degree=2)
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
