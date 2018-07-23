from dolfin import *
from FEMModel import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2017/11/16"


class UncertainBoundariesDiffusionFEMML(FEMModel):
    def __init__(self, a, f, bc_var, M_gen, mesh_size):
        self.a         = a
        self.f         = f
        self.bc_var    = bc_var
        self.M_gen     = M_gen

        self.mesh_size = mesh_size
        self.init_simple_mesh()

    def solve(self, z):
        # Make FEniCS output only the most important messages
        set_log_level(WARNING)

        # Create mesh if there is none
        if not hasattr(self, 'mesh'):
            self.init_simple_mesh()

        # Create approximation space
        V = FunctionSpace(self.mesh, 'Lagrange', 2)

        params = self.split_params([self.a, self.f, self.bc_var], z)

        # Define boundary conditions
        
	tol = 1E-14
	u_zero = Constant(0.0) # Defines a constant boundary condition for the top and bottom parts of the rectangle
	def boundary_T(x, on_boundary): # Boundary top
	    return on_boundary and near(x[1],0,tol)
	
	def boundary_B(x, on_boundary): # Boundary bottom
	    return on_boundary and near(x[1],1,tol)
	bc_T = DirichletBC(V,u_zero, boundary_T)
	bc_B = DirichletBC(V,u_zero, boundary_B)
	
	# Define 4 (why 4? Because that's a fair amount ... completely random) piecewise constant left boundary conditions
	separator = 1.0/4.0
	def boundary_L1(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],0,tol) and x[1] >= 0-tol and x[1] <= separator + tol
	bc_L1 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][0]), boundary_L1)
	def boundary_L2(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],0,tol) and x[1] >= separator-tol and x[1] <= 2*separator + tol
	bc_L2 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][1]), boundary_L2)
	def boundary_L3(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],0,tol) and x[1] >= 2*separator-tol and x[1] <= 3*separator + tol
	bc_L3 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][2]), boundary_L3)
	def boundary_L4(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],0,tol) and x[1] >= 3*separator-tol and x[1] <= 1 + tol
	bc_L4 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][3]), boundary_L4)
	# And now 4 BC for the right boundary
	def boundary_R1(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],1,tol) and x[1] >= 0-tol and x[1] <= separator + tol
	bc_R1 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][4]), boundary_R1)
	def boundary_R2(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],1,tol) and x[1] >= separator-tol and x[1] <= 2*separator + tol
	bc_R2 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][5]), boundary_R2)
	def boundary_R3(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],1,tol) and x[1] >= 2*separator-tol and x[1] <= 3*separator + tol
	bc_R3 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][6]), boundary_R3)
	def boundary_R4(x,on_boundary): # Boundary top
	    return on_boundary and near(x[0],1,tol) and x[1] >= 3*separator-tol and x[1] <= 1 + tol
	bc_R4 = DirichletBC(V, Constant(self.bc_var.c+self.bc_var.local_var*params[2][7]), boundary_R4)
	bcs = [bc_T, bc_B, bc_L1, bc_L2, bc_L3, bc_L4, bc_R1, bc_R2, bc_R3, bc_R4]

        # Define variational problem
        w = TrialFunction(V)
        v = TestFunction(V)


        x = SpatialCoordinate(self.mesh)
        A = self.a(x, Constant(params[0])) * inner(nabla_grad(w), nabla_grad(v)) * dx
        L = self.f(x, Constant(params[1])) * v * dx
        # Create goal-functional for error estimation
        u      = Function(V)
        self.M = self.M_gen(self, u, dx)

        # Create solver
        problem     = LinearVariationalProblem(A, L, u, bcs)
	self.solver = LinearVariationalSolver(problem) #, solver_parameters={'linear_solver': 'iterative'})
        #self.solver.parameters["linear_solver"] ="iterative"
	# y[k] = assemble(myAverage(mesh, u, dx))

        # Compute solution
        self.solver.solve()

        return u
