include("../ArcticOcean.jl")

using .ArcticOcean
using .ArcticOcean: grid, ocean, atmosphere
using .ArcticOcean: SI_meta_init, SC_meta_init
using .ArcticOcean: arctic_outputs!

#####
##### A Prognostic Sea-ice model
#####

using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

# Remember to pass the SSS as a bottom bc to the sea ice!
SSS = view(ocean.model.tracers.S.data, :, :, grid.Nz)
bottom_heat_boundary_condition = IceWaterThermalEquilibrium(SSS)

SSU = view(ocean.model.velocities.u, :, :, grid.Nz)
SSV = view(ocean.model.velocities.u, :, :, grid.Nz)
τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV)
τua = Field{Face, Center, Nothing}(grid)
τva = Field{Center, Face, Nothing}(grid)

dynamics = SeaIceMomentumEquation(grid;
                                  coriolis = ocean.model.coriolis,
                                  top_momentum_stress = (u=τua, v=τva),
                                  bottom_momentum_stress = τo,
                                  ocean_velocities = (u=0.1*SSU, v=0.1*SSV),
                                  rheology = ElastoViscoPlasticRheology(),
                                  solver = SplitExplicitSolver(120))

sea_ice = ClimaOcean.SeaIceSimulations.sea_ice_simulation(grid; bottom_heat_boundary_condition, dynamics, advection=WENO(order=7))
set!(sea_ice.model.ice_thickness,     SI_meta_init; inpainting=nothing)
set!(sea_ice.model.ice_concentration, SC_meta_init; inpainting=nothing)

arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere)
arctic = Simulation(arctic, Δt=2minutes, stop_time=365days)

ArcticOcean.arctic_outputs!(arctic, "EVP-rheology/")

# And add it as a callback to the simulation.
add_callback!(arctic, ArcticOcean.progress, IterationInterval(10))

run!(arctic)
