module ArcticOcean

using ClimaOcean
using PythonCall
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Dates

const start_date = DateTime(1996, 11, 1)
const end_date   = DateTime(1997, 06, 1)

const arch = GPU()
const Nx = 3600
const Ny = 400
const Nz = 50

include("arctic_grid.jl")
include("glorys_restoring.jl")
include("atmosphere.jl")
include("ocean.jl")
include("outputs.jl")

end
