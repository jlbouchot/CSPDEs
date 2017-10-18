from fenics import *

# Create the mesh and the functions spaces on it
nx = 50
ny = nx
nz = nx
myMesh = UnitCubeMesh(nx, ny, nz) # defined via the number of points in each directions
V = FunctionSpace(myMesh, 'Lagrange', 1) # Degree 1 Lagrange interpolating polynomials

# Define the boundary conditions
# u_boundary = Constant(10)
u_boundary = Expression('1 + x[0]*x[0] + 2*x[1]*x[1] + 3*x[2]*x[2]', degree=2) # Make sure we have a good enough approximation of this function. 

def boundary(x, on_boundary): 
    return on_boundary # a predefined boolean function taking the mathematical boundaries of the domain

## Alternate approaches
# tol_boundary = 1e-10
# def boundary(x):
#     return abs(x[0]) < tol_boundary or abs(x[1]) < tol_boundary or abs(x[2]) < tol_boundary or abs(x[0]-1) < tol_boundary or abs(x[1]-1) < tol_boundary or abs(x[2]-1) < tol_boundary
# # other option: use the near(x,1, tolerance) function
bc = DirichletBC(V, u_boundary, boundary)

# Set up the variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-12.0)
# The bilinear form
a = dot(grad(u), grad(v))*dx
# The right hand side linear form
L = f*v*dx

# Et voila! 
u = Function(V)
solve(a == L, u, bc)


plot(u)
plot(myMesh)
interactive()
