DO_DOLFINX = False

if DO_DOLFINX:
    import dolfinx
    import ufl
    from mpi4py import MPI
    from dolfinx.fem.petsc import LinearProblem
    import pyvista

else:
    from dolfin import *
    import matplotlib.pyplot as plt

import sys
import numpy as np
import argparse

import time
import resource

import pandas as pd


__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2025, INRIA, RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich, and Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.5.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2025/07/01"
__creationdate__ = "2025/06/17"




# def Main(outfile, d = 10, L_max = 4, orig_mesh_size = 2000):
def Main(outfile = "thatTest", d = 5, grid_points = tuple([2000, 2000]), L_max = 4, algo_name = "whtp", gamma = 1.035, L_min = 1, sampling_name = "p", nb_iter = 500, epsilon = 1e-3, nb_tests = None, mu = 2.0, abar = 4.3, imp = 1, w_cst = 0.5, dat_constant = 10, experiment_name = "weighted_cosine_avg_p_2D", tensor_based=True, ansatz_space = 0, t_0 = 1, t_prime = 1, p0 = 1./4., p = 3./10., const_sJ = 5, no_compute=False, exponent=1.0/4.0):


    dict_config = {'d': d, 'J': L_min, "L": L_max, "h0": grid_points, "vj": gamma, "weightCosine":w_cst, "nbSamples": sampling_name, "Tensor": tensor_based, 't': t_0, "tprime": t_prime, 'p0': p0, "p": p, "s_J": const_sJ, "s_L": dat_constant, "trig_power": mu, "abar":abar, "energy_fluctuations": imp, "algo": algo_name, "iter": nb_iter, "tolres": epsilon, "ansatz": ansatz_space, "no_compute": no_compute, "alpha": exponent}


    # Create mesh if there is none
    h0 = 20


    nb_tests = 1
    nb_levels = 2
    L_start = 0
    preconditioners = {"gmres": ["amg", "ilu"]}
    # preconditioners = {"gmres": ["amg", "petsc_amg", "hypre_amg"], "petsc": ["lu"]}
    # preconditioners = {"gmres": ["amg", "petsc_amg", "hypre_amg"], "petsc": ["lu", "SuperLU"]}
    preonly = {"petsc": True}

    results = pd.DataFrame(columns=["solver", "preconditioner", "level", "wall-clock", "user", "cpu"])

    # for s in ["petsc", "gmres"]:
    all_wc_time = time.time()
    all_s_resource = resource.getrusage(resource.RUSAGE_SELF)
    # for s in ["mumps", "superlu", "petsc", "umfpack", "gmres"]:
    for s in ["petsc"]: 
        print(f"Current solver is {s}")

        if s in preconditioners.keys():
            pcs = preconditioners[s]
        else:
            pcs = [None]

        if s in preonly.keys():
            ksp = preonly[s]
        else:
            ksp = False
        
        for p in pcs:
            print(f"Current preconditioner is {p}")
            for j in range(L_start,nb_levels):
                user_time = 0
                wc_time = 0
                cpu_time = 0

                if not DO_DOLFINX:
                    for i in range(nb_tests): 
                        # FENICS 
                        # Make FEniCS output only the most important messages
                        CRITICAL  = 50 #, // errors that may lead to data corruption and suchlike
                        ERROR     = 40 #, // things that go boom
                        WARNING   = 30 #, // things that may go boom later
                        INFO      = 20 #, // information of general interest
                        PROGRESS  = 16 #, // what's happening (broadly)
                        TRACE     = 13 #, // what's happening (in detail)
                        DBG       = 10#   // sundry
                        set_log_level(ERROR)
                        mesh = UnitSquareMesh(h0*2**j,h0*2**j)

                        # Create approximation space
                        V = FunctionSpace(mesh, 'Lagrange', 1)
                        # Define boundary conditions
                        bc = DirichletBC(V, Constant(0.0), lambda x, on_boundary: on_boundary)
                        # Define variational problem
                        w = TrialFunction(V)
                        v = TestFunction(V)
                        cst1 = 0.5
                        cst2 = 0.25
                        x = SpatialCoordinate(mesh)
                        A = inner(nabla_grad(w), nabla_grad(v)) * dx
                        L = v * dx
                        # Create goal-functional for error estimation
                        u      = Function(V)
                        # Create solver
                        problem     = LinearVariationalProblem(A, L, u, bc)
                        solver = LinearVariationalSolver(problem) 

                    

                        ############################################
                        ## Great candidate so far: 
                        # solver.parameters["linear_solver"] = "gmres"
                        # solver.parameters["preconditioner"] = "amg"
                        ############################################
                        ############################################
                        # Good for wall clock time, not so much for CPU
                        # solver.parameters["linear_solver"] = "mumps"
                        ############################################
                        solver.parameters["linear_solver"] = s
                        if p is not None: 
                            solver.parameters["preconditioner"] = p

                        # if ksp: 
                        #     solver.parameters["ksp_type"] = "preonly"
                        
                        # solver.parameters["linear_solver"] = "iterative"
                        # solver.parameters["preconditioner"] = "amg"
                        # solver.parameters["linear_solver"] = "superlu_dist"
                        # solver.parameters["linear_solver"] = "superlu"
                        # solver.parameters["preconditioner"] = "amg" 
                        # solver.parameters["preconditioner"] = "amg"
                        #solver.parameters["preconditioner"] = "amg"
                        # solver.parameters["preconditioner"] = "ilu"
                        # solver.parameters["linear_solver"] = "petsc"
                        # solver.parameters["linear_solver"] = "iterative" -> This is garbage without further parameters
                        # solver.parameters["preconditioner"] = "petsc_amg"
                        # solver.parameters["relative_tolerance"] = 1e-3
                        # solver.parameters["absolute_tolerance"] =1e-6
                        # solver.parameters["linear_solver"] ="iterative"
                        # solver.parameters["linear_solver"] = "umfpack"
                        # solver.parameters["preconditioner"] = "petsc_amg"
                        # y[k] = assemble(myAverage(mesh, u, dx))
                        # Compute solution
                        tWC = time.time()
                        start = resource.getrusage(resource.RUSAGE_SELF)
                        solver.solve()
                        end = resource.getrusage(resource.RUSAGE_SELF)
                        wc_time += time.time() - tWC 
                        u_time = end.ru_utime - start.ru_utime
                        user_time += u_time
                        c_time = end.ru_stime - start.ru_stime + user_time
                        cpu_time += c_time



                else: 
                    # FENICSX
                    # mesh = dolfinx.mesh.create_unit_square(h0*2**j,h0*2**j)
                    for i in range(nb_tests): 
                        mesh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, h0*2**j,h0*2**j, dolfinx.mesh.CellType.triangle)

                        # Create approximation space
                        V = dolfinx.fem.functionspace(mesh, ('Lagrange', 1))
                        # Define boundary conditions
                        boundary_dofs = dolfinx.fem.locate_dofs_geometrical(V, on_boundary_dolfinx)
                        bc = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs, V)
                        # bc = dolfinx.fem.dirichletbc(V, dolfinx.fem.Constant(mesh, dolfinx.default_scalar_type(0.0)), lambda x, on_boundary: on_boundary)
                        # Define variational problem
                        w = ufl.TrialFunction(V)
                        v = ufl.TestFunction(V)
                        cst1 = 0.5
                        cst2 = 0.25
                        # x = SpatialCoordinate(mesh)
                        A = ufl.dot(ufl.grad(w), ufl.grad(v)) * ufl.dx
                        L = v * ufl.dx
                        # Create solver
                        if p == "lu":
                            problem     = dolfinx.fem.petsc.LinearProblem(A, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": p})
                        else:
                            problem     = dolfinx.fem.petsc.LinearProblem(A, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "superlu"})
                        # solver = LinearVariationalSolver(problem)
                        tWC = time.time()
                        start = resource.getrusage(resource.RUSAGE_SELF)
                        problem.solve()
                        end = resource.getrusage(resource.RUSAGE_SELF)
                        wc_time += time.time() - tWC 
                        u_time = end.ru_utime - start.ru_utime
                        user_time += u_time
                        c_time = end.ru_stime - start.ru_stime + user_time
                        cpu_time += c_time 

                

                results = pd.concat([pd.DataFrame([[s, p, j, wc_time, user_time, cpu_time]], columns=results.columns), results], ignore_index=True)

        # print(f" solver {s} has parameters {dict(solver.parameters)}")

                    

    print(results)
    results.to_csv("BenchmarkingSolvers.csv")
    all_wc_time = time.time() - all_wc_time
    all_e_resource = resource.getrusage(resource.RUSAGE_SELF)
    total_u_time = all_e_resource.ru_utime - all_s_resource.ru_utime
    total_c_time = all_e_resource.ru_stime - all_s_resource.ru_stime + total_u_time
    print(f"All experiments ran in WCTime: {all_wc_time}s -- CPUTime: {total_c_time}s -- USERTime: {total_u_time}")
    if not DO_DOLFINX:
        list_linear_solver_methods()
        list_krylov_solver_preconditioners()

    # pyvista.start_xvfb(1.0)
    plot_mesh(mesh, L = nb_levels-1)
    

def on_boundary_dolfinx(x):
    return np.logical_or(np.logical_or(np.isclose(x[0], 0), np.isclose(x[0], 1)),np.logical_or(np.isclose(x[1], 0), np.isclose(x[1], 1)))


def plot_mesh(mesh, values = None, L = None):
# def plot_mesh(mesh: dolfinx.mesh.Mesh, values = None):
    """
    Given a DOLFINx mesh, create a `pyvista.UnstructuredGrid`,
    and plot it and the mesh nodes.

    Args:
        mesh: The mesh we want to visualize
        values: List of values indicating a marker for each cell in the mesh

    Note:
        If `values` are given as input, they are assumed to be a marker
        for each cell in the domain.
    """
    # We create a pyvista plotter instance
    if DO_DOLFINX:
        plotter = pyvista.Plotter()

        # # Since the meshes might be created with higher order elements,
        # # we start by creating a linearized mesh for nicely inspecting the triangulation.
        # V_linear = dolfinx.fem.functionspace(mesh, ('Lagrange', 1))
        # # V_linear = FunctionSpace(mesh, "Lagrange", 1)
        # linear_grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V_linear))

        # If the mesh is higher order, we plot the nodes on the exterior boundaries,
        # as well as the mesh itself (with filled in cell markers)
        if mesh.geometry.cmap.degree > 1:
            ugrid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(mesh))
            if values is not None:
                ugrid.cell_data["Marker"] = values
            plotter.add_mesh(ugrid, style="points", color="b", point_size=10)
            ugrid = ugrid.tessellate()
            plotter.add_mesh(ugrid, show_edges=False)
            plotter.add_mesh(linear_grid,style="wireframe", color="black")
        else:
            # If the mesh is linear we add in the cell markers
            if values is not None:
                linear_grid.cell_data["Marker"] = values
            grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(mesh))
            plotter.add_mesh(grid,show_edges=True)

        # We plot the coordinate axis and align it with the xy-plane
        plotter.show_axes()
        plotter.view_xy()
        # plotter.show(screenshot="mesh_plot.png")
        if not pyvista.OFF_SCREEN:
            plotter.show()  
    else: 
        # Plot with matplotlib
        plot(mesh)
        fig_title =  f"Triangular mesh L ={L}" if L else f"Triangular mesh"
        plt.title(fig_title)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.show()

# import matplotlib.pyplot as plt
# from dolfinx.plot import create_vtk_topology, extract_geometry

# topology, cell_types = create_vtk_topology(mesh, mesh.topology.dim)
# points = extract_geometry(mesh) 

### Main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-o", "--output-file", help="File to write the results", default="outputDiffusionMLPolynomial", required=False)
    parser.add_argument("-L", "--nb-level", help="Number of levels used", default=4, required=False)
    parser.add_argument("-x", "--mesh-x", help="Size of the coarsest level (number of grid points) in the x direction", default=2000, required=False)
    parser.add_argument("-y", "--mesh-y", help="Size of the coarsest level (number of grid points) in the y direction", default=2000, required=False)
    parser.add_argument("-N", "--nb-iter", help="Number of iterations for the (potential) iterative greedy algorithm", default=50, required=False)
    parser.add_argument("-r", "--recovery-algo", help="String for the algorithm for weighted l1 recovery", default="whtp", required=False)
    parser.add_argument("-s", "--l-start", help="Instead of going through all the levels, give it a starting point", default=1, required=False)
    parser.add_argument("-t", "--sampling", help="Select a sampling strategy (pragmatic or theoretic or new)", default="pragmatic", required=False)
    parser.add_argument("-n", "--nb-tests", help="Number of tests 'on the fly'", default=None, required=False)
    parser.add_argument("-p", "--power", help="Power of the decay of the trigonometric expansion (~ mu)", default=4.0, required=False)
    parser.add_argument("-a", "--abar", help="Value of the mean field", default=10, required=False)
    parser.add_argument("-b", "--better-compute", help="Should the computations be done on the fly, using tensor representation (Default is TRUE)", default="True", required=False)
    parser.add_argument("-c", "--dat_constant", help="Multiplicative constant for expression of s_L", default=15., required=False)
    parser.add_argument("-d", "--nb-cosines", help="Number of random cosine and sine parameters", default=5, required=False)
    parser.add_argument("-e", "--tol-res", help="Tolerance on the residual for the recovery algorithms (called epsilon everywhere)", default=1e-4, required=False)
    parser.add_argument("-E", "--exponent", help="Power of the polynomial weight (~ alpha)", default=1.0/4.0, required=False)
    parser.add_argument("-f", "--prefix-precompute", help="How should the precomputed data for this test be called?", default="testingWCosine", required=False)
    parser.add_argument("-g", "--gamma", help="Value of the constant coefficients", default=1.035, required=False)
    parser.add_argument("-i", "--fluctuation-importance", help="What is the importance of the fluctuations with respect to the mean field (default is 1)", default=1, required=False)
    parser.add_argument("-j", "--ansatz-space", help="What type of Ansatz space is used? (Default is 0)", default="0", required=False)
    parser.add_argument("-k", "--no-compute", help="Should we skip all computations and only check values for s, m, and N (Default = False)", default=False, required=False)
    parser.add_argument("-w", "--weight-cosine", help="How much weight the local cosine carries (Default = 1, ~ Upsilon)", default=1, required=False)
    parser.add_argument("--t_0", help="What is the smoothness of the data (Default is 1)", default="1", required=False)
    parser.add_argument("--t_prime", help="What is the smoothness of the functional (Default is 1)", default="1", required=False)
    parser.add_argument("--smooth_0", help="What kind of smoothness in the original space can be expected (Default is 1/2)", default="0.5", required=False)
    parser.add_argument("--smooth_t", help="What kind of smoothness in the smooth space can be expected (Default is 1/2)", default="0.5", required=False)
    parser.add_argument("--const_sJ", help="What is the expected constant in the expression of s_J (Default is 15)", default="15", required=False)
    args = parser.parse_args()
	
    
    Main(args.output_file, int(args.nb_cosines), tuple([int(args.mesh_x),int(args.mesh_y)]), int(args.nb_level), args.recovery_algo.lower(), float(args.gamma), int(args.l_start), args.sampling, int(args.nb_iter), float(args.tol_res), None if args.nb_tests is None else int(args.nb_tests), float(args.power), float(args.abar), float(args.fluctuation_importance), float(args.weight_cosine), float(args.dat_constant), args.prefix_precompute, args.better_compute.lower()=="true", int(args.ansatz_space), float(args.t_0), float(args.t_prime), float(args.smooth_0), float(args.smooth_t), float(args.const_sJ), False if args.nb_tests is None else args.no_compute, float(args.exponent))
    # Main(sys.argv[1])

