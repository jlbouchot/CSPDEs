from dolfinx import mesh, plot
from mpi4py import MPI
import pyvista

domain = mesh.create_unit_square(MPI.COMM_WORLD, 4, 4, mesh.CellType.triangle)
# topology, cell_types, geometry = plot.vtk_mesh(domain)
# grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(domain))
plotter = pyvista.Plotter()
plotter.add_mesh(grid, show_edges=True)
plotter.show()
