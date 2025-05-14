module ArcticOcean

using ClimaOcean
using PythonCall
using ClimaSeaIce
using Oceananigans
using Dates

const start_date = Date(1996, 11, 1)
const end_date   = Date(1997, 06, 1)

const arch = GPU()
const Nx = 4320
const Ny = 450
const Nz = 60

include("arctic_grid.jl")
include("glorys_restoring.jl")
include("atmosphere.jl")
include("ocean.jl")
include("outputs.jl")

end