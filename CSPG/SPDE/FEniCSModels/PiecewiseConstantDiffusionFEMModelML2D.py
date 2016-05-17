from dolfin import *
from FEMModel import *
import numpy as np

# Still have to find the way to 
# 1) Automatize the process to an unknown d number of sub domains and 
# 2) uncertain splitting points
class Omega00(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 0 and x[0] <= 1.0/5.0 and x[1] >= 0 and x[1] <= 1.0/5.0 else False

class Omega01(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 0 and x[0] <= 1.0/5.0 and x[1] >= 1.0/5.0 and x[1] <= 2.0/5.0 else False

class Omega02(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 0 and x[0] <= 1.0/5.0 and x[1] >= 2.0/5.0 and x[1] <= 3.0/5.0 else False
		
class Omega03(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 0 and x[0] <= 1.0/5.0 and x[1] >= 3.0/5.0 and x[1] <= 4.0/5.0 else False
		
class Omega04(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 0 and x[0] <= 1.0/5.0 and x[1] >= 4.0/5.0 and x[1] <= 1.0 else False

class Omega10(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 1.0/5.0 and x[0] <= 2.0/5.0 and x[1] >= 0 and x[1] <= 1.0/5.0 else False

class Omega11(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 1.0/5.0 and x[0] <= 2.0/5.0 and x[1] >= 1.0/5.0 and x[1] <= 2.0/5.0 else False

class Omega12(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 1.0/5.0 and x[0] <= 2.0/5.0 and x[1] >= 2.0/5.0 and x[1] <= 3.0/5.0 else False
		
class Omega13(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 1.0/5.0 and x[0] <= 2.0/5.0 and x[1] >= 3.0/5.0 and x[1] <= 4.0/5.0 else False
		
class Omega14(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 1.0/5.0 and x[0] <= 2.0/5.0 and x[1] >= 4.0/5.0 and x[1] <= 1.0 else False
		
class Omega20(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 2.0/5.0 and x[0] <= 3.0/5.0 and x[1] >= 0 and x[1] <= 1.0/5.0 else False

class Omega21(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 2.0/5.0 and x[0] <= 3.0/5.0 and x[1] >= 1.0/5.0 and x[1] <= 2.0/5.0 else False

class Omega22(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 2.0/5.0 and x[0] <= 3.0/5.0 and x[1] >= 2.0/5.0 and x[1] <= 3.0/5.0 else False
		
class Omega23(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 2.0/5.0 and x[0] <= 3.0/5.0 and x[1] >= 3.0/5.0 and x[1] <= 4.0/5.0 else False
		
class Omega24(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 2.0/5.0 and x[0] <= 3.0/5.0 and x[1] >= 4.0/5.0 and x[1] <= 1.0 else False
		
class Omega30(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 3.0/5.0 and x[0] <= 4.0/5.0 and x[1] >= 0 and x[1] <= 1.0/5.0 else False

class Omega31(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 3.0/5.0 and x[0] <= 4.0/5.0 and x[1] >= 1.0/5.0 and x[1] <= 2.0/5.0 else False

class Omega32(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 3.0/5.0 and x[0] <= 4.0/5.0 and x[1] >= 2.0/5.0 and x[1] <= 3.0/5.0 else False
		
class Omega33(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 3.0/5.0 and x[0] <= 4.0/5.0 and x[1] >= 3.0/5.0 and x[1] <= 4.0/5.0 else False
		
class Omega34(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 3.0/5.0 and x[0] <= 4.0/5.0 and x[1] >= 4.0/5.0 and x[1] <= 1.0 else False
		
class Omega40(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 4.0/5.0 and x[0] <= 1.0 and x[1] >= 0 and x[1] <= 1.0/5.0 else False

class Omega41(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 4.0/5.0 and x[0] <= 1.0 and x[1] >= 1.0/5.0 and x[1] <= 2.0/5.0 else False

class Omega42(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 4.0/5.0 and x[0] <= 1.0 and x[1] >= 2.0/5.0 and x[1] <= 3.0/5.0 else False
		
class Omega43(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 4.0/5.0 and x[0] <= 1.0 and x[1] >= 3.0/5.0 and x[1] <= 4.0/5.0 else False
		
class Omega44(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] >= 4.0/5.0 and x[0] <= 1.0 and x[1] >= 4.0/5.0 and x[1] <= 1.0 else False

tol = 1E-14   # tolerance for coordinate comparisons
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)
		
# class BottomBoundary(SubDomain):
    # def inside(self, x, on_boundary):
        # return on_boundary and abs(x) < tol


# class TopBoundary(SubDomain):
    # def inside(self, x, on_boundary):
        # return on_boundary and abs(x - 1) < tol


class PiecewiseConstantDiffusionFEMModelML2D(FEMModel):
    
    def __init__(self, abar, f, var, mesh_size, M_gen):
        
        self.abar      = abar
        self.f         = f
        self.var       = var
        self.M_gen     = M_gen
        self.num_subd  = 24

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
        subdomains = MeshFunction('size_t', self.mesh, 2) # The last argument corresponds to the dimension of the cells: here, intervals, dim 1. 
        subdomains.set_all(0) ### IMPORTANT! This is not a GOOD way to deal with it. But I have no clue how to solve the problem. It appears that the 'meshes' at the boundaries between two partitions are marked with random numbers
        subdomain0 = Omega00()
        subdomain0.mark(subdomains, 0)
        subdomain1 = Omega01()
        subdomain1.mark(subdomains, 1)
        subdomain2 = Omega02()
        subdomain2.mark(subdomains, 2)
        subdomain3 = Omega03()
        subdomain3.mark(subdomains, 3)
        subdomain4 = Omega04()
        subdomain4.mark(subdomains, 4)
        subdomain0 = Omega00()
        subdomain0.mark(subdomains, 5)
        subdomain1 = Omega01()
        subdomain1.mark(subdomains, 6)
        subdomain2 = Omega02()
        subdomain2.mark(subdomains, 7)
        subdomain3 = Omega03()
        subdomain3.mark(subdomains, 8)
        subdomain4 = Omega04()
        subdomain4.mark(subdomains, 9)
		subdomain5 = Omega10()
        subdomain5.mark(subdomains, 10)
        subdomain6 = Omega11()
        subdomain6.mark(subdomains, 11)
        subdomain7 = Omega12()
        subdomain7.mark(subdomains, 12)
        subdomain8 = Omega13()
        subdomain8.mark(subdomains, 13)
        subdomain9 = Omega14()
        subdomain9.mark(subdomains, 14)
		subdomain10 = Omega20()
        subdomain10.mark(subdomains, 10)
        subdomain11 = Omega21()
        subdomain11.mark(subdomains, 11)
        subdomain12 = Omega22()
        subdomain12.mark(subdomains, 12)
        subdomain13 = Omega23()
        subdomain13.mark(subdomains, 13)
        subdomain14 = Omega24()
        subdomain14.mark(subdomains, 14)
		subdomain15 = Omega30()
        subdomain15.mark(subdomains, 15)
        subdomain16 = Omega31()
        subdomain16.mark(subdomains, 16)
        subdomain17 = Omega32()
        subdomain17.mark(subdomains, 17)
        subdomain18 = Omega33()
        subdomain18.mark(subdomains, 18)
        subdomain19 = Omega34()
        subdomain19.mark(subdomains, 19)
		subdomain20 = Omega40()
        subdomain20.mark(subdomains, 20)
        subdomain21 = Omega41()
        subdomain21.mark(subdomains, 21)
        subdomain22 = Omega42()
        subdomain22.mark(subdomains, 22)
        subdomain23 = Omega43()
        subdomain23.mark(subdomains, 23)
        subdomain24 = Omega44()
        subdomain24.mark(subdomains, 24)

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
        V = FunctionSpace(self.mesh, 'Lagrange', degree=2)
        # Gamma_0 = DirichletBC(V, Constant(0), BottomBoundary())
        # Gamma_1 = DirichletBC(V, Constant(1), TopBoundary())
		Gamma_l = DirichletBC(V, Constant(0), Left())
        Gamma_r = DirichletBC(V, Constant(1), Right())
		Gamma_t = DirichletBC(V, Constant(0), Top())
        Gamma_b = DirichletBC(V, Constant(1), Bottom())

        # Define boundary conditions
        bc = [Gamma_l, Gamma_r, Gamma_b, Gamma_t] # DirichletBC(V, Constant(0.0), lambda x, on_boundary: on_boundary)

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