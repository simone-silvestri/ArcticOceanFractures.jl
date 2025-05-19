using Oceananigans
using ClimaOcean

r_faces = ClimaOcean.exponential_z_faces(; Nz, depth=6000)
z_faces = MutableVerticalDiscretization(r_faces)

grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces, southernmost_latitude=40, halo=(7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=20, major_basins=10, interpolation_passes=10)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
