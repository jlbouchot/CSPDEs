import sys, numpy, scipy, numba, cvxpy, matplotlib, pandas, progressbar, mpi4py, petsc4py
from dolfin import *
print("python", sys.version.split()[0], "| dolfin", __import__("dolfin").__version__, "| numpy", numpy.__version__, "| scipy", scipy.__version__, "| cvxpy", cvxpy.__version__)
mesh = UnitSquareMesh(32, 32); V = FunctionSpace(mesh, "Lagrange", 1)
u, v = TrialFunction(V), TestFunction(V)
a = inner(Constant(10.0)*grad(u), grad(v))*dx; L = Constant(1.0)*v*dx
bc = DirichletBC(V, Constant(0.0), "on_boundary")
for solver in [("petsc", "ilu"), ("cg", "amg"), ("gmres", "ilu"), ("superlu_dist", "none"), ("umfpack", "none"), ("mumps", "none")]:
    w = Function(V); solve(a == L, w, bc, solver_parameters={"linear_solver": solver[0], "preconditioner": solver[1]} if solver[1] != "none" else {"linear_solver": solver[0]})
    print("solve", solver, "avg =", assemble(w*dx))
x = cvxpy.Variable(3); p = cvxpy.Problem(cvxpy.Minimize(cvxpy.norm1(x)), [x[0] + x[1] == 1]); p.solve(); print("cvxpy", p.status, round(p.value, 6))
print("MPI size", mpi4py.MPI.COMM_WORLD.Get_size())
