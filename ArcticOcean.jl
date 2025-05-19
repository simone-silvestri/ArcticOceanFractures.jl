module ArcticOcean

using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Dates

const start_date  = DateTime(1996, 1,   1)
const end_date    = DateTime(1997, 12, 31)

const arch = GPU()
const Nx = 3600
const Ny = 500
const Nz = 50

include("arctic_grid.jl")
include("ecco_restoring.jl")
include("atmosphere.jl")
include("ocean.jl")
include("outputs.jl")

end
