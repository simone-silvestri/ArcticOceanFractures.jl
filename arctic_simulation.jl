using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf

using CUDA

start_date = Date(1996, 11, 1)
end_date   = Date(1998, 11, 1)

##### A Prognostic Ocean model
#####

#####
##### A Prognostic Sea-ice model
#####

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

sea_ice = sea_ice_simulation(grid; bottom_heat_boundary_condition, dynamics, advection=WENO(order=7))

#####
##### A Prescribed Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(40), dataset=MultiyearJRA55())
radiation  = Radiation()

#####
##### Arctic coupled model
#####

arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
arctic = Simulation(arctic, Δt=2minutes, stop_time=365days)

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    hmax  = maximum(sea_ice.model.ice_thickness)
    ℵmax  = maximum(sea_ice.model.ice_concentration)
    uimax = maximum(abs, sea_ice.model.velocities.u)
    vimax = maximum(abs, sea_ice.model.velocities.v)
    uomax = maximum(abs, ocean.model.velocities.u)
    vomax = maximum(abs, ocean.model.velocities.v)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg3 = @sprintf("max uᵢ: (%.2f, %.2f) m s⁻¹, ", uimax, vimax)
    msg4 = @sprintf("max uₒ: (%.2f, %.2f) m s⁻¹, ", uomax, vomax)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(arctic, progress, IterationInterval(10))

run!(arctic)

#####
##### Comparison to ECCO Climatology
#####

dataset = ECCO4Monthly()
dates   = all_dates(version)[1:12]

h_metadata = Metadata(:sea_ice_thickness;     dataset, dates)
ℵ_metadata = Metadata(:sea_ice_concentration; dataset, dates)

# Montly averaged ECCO data
hE = ECCOFieldTimeSeries(h_metadata, grid; time_indices_in_memory=12)
ℵE = ECCOFieldTimeSeries(ℵ_metadata, grid; time_indices_in_memory=12)

# Daily averaged Model output
h = FieldTimeSeries("averaged_sea_ice_quantities.jld2", "h")
ℵ = FieldTimeSeries("averaged_sea_ice_quantities.jld2", "ℵ")

# Montly average the model output
hm = FieldTimeSeries{Center, Center, Nothing}(grid, hE.times; backend=InMemory())
ℵm = FieldTimeSeries{Center, Center, Nothing}(grid, hE.times; backend=InMemory())

hMe = [mean(h[t]) for t in 1:length(h.times)]
ℵMe = [mean(ℵ[t]) for t in 1:length(ℵ.times)]

