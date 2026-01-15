import numpy as np
import petsc4py
import petsc4py.PETSc

import ufl
import dolfinx 
import dolfinx.fem
import dolfinx.fem.petsc
import dolfinx.mesh
import dolfinx.io

import mpi4py
import mpi4py.MPI

import time
import resource


h = 1280 # Corresponds to L = 6, with h0 = 20
mesh = dolfinx.mesh.create_unit_square(mpi4py.MPI.COMM_WORLD, h,h, dolfinx.mesh.CellType.triangle)

V = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

dx = ufl.dx

f = dolfinx.fem.Function(V)
dim_V = f.x.index_map.size_local
# assert dim_V == 4
f.x.array[:] = np.arange(1, dim_V + 1)

a = ufl.inner(u, v) * dx
F = ufl.inner(f, v) * dx
a_cpp = dolfinx.fem.form(a)
F_cpp = dolfinx.fem.form(F)

A = dolfinx.fem.petsc.assemble_matrix(a_cpp)
b = dolfinx.fem.petsc.assemble_vector(F_cpp)

A.assemble()

solution = dolfinx.fem.Function(V)

for (package_number, package) in enumerate(("umfpack", None, "mumps", "superlu", "superlu_dist")):
    ksp = petsc4py.PETSc.KSP().create()
    tWC = time.time()
    start = resource.getrusage(resource.RUSAGE_SELF)
    ksp.setOperators(A)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    if package is not None:
        ksp.getPC().setFactorSolverType(package)
    ksp.setFromOptions()
    ksp.solve(b, solution.x.petsc_vec)
    solution.x.petsc_vec.ghostUpdate(addv=petsc4py.PETSc.InsertMode.INSERT, mode=petsc4py.PETSc.ScatterMode.FORWARD)
    assert np.allclose(solution.x.array, np.arange(1, dim_V + 1))
    # with dolfinx.io.VTXWriter(mesh.comm, "solution.bp", solution) as vtx_file:
    #     vtx_file.write(package_number * 1.0)
    end = resource.getrusage(resource.RUSAGE_SELF)
    wc_time = time.time() - tWC 
    u_time = end.ru_utime - start.ru_utime
    print(f"Wall clock time for solver {package} is {wc_time}")
    print(f"User time for solver {package} is {u_time}")
    c_time = end.ru_stime - start.ru_stime + u_time
    print(f"CPU time for solver {package} is {c_time}")
    ksp.destroy()