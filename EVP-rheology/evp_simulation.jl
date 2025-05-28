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
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Operators: Δx⁻¹ᶠᶜᶜ, Δy⁻¹ᶜᶠᶜ
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies
using ClimaSeaIce: IceWaterThermalEquilibrium

SSS = interior(ocean.model.tracers.S, :, :, grid.Nz:grid.Nz)
bottom_heat_boundary_condition = IceWaterThermalEquilibrium(SSS)

SSU = view(ocean.model.velocities.u, :, :, grid.Nz)
SSV = view(ocean.model.velocities.v, :, :, grid.Nz)
coriolis = ocean.model.coriolis

τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV, Cᴰ=sea_ice_ocean_drag_coefficient)
τua = Field{Face, Center, Nothing}(grid)
τva = Field{Center, Face, Nothing}(grid)

dynamics = SeaIceMomentumEquation(grid;
                                  coriolis,
                                  top_momentum_stress = (u=τua, v=τva),
                                  bottom_momentum_stress = τo,
                                  solver = SplitExplicitSolver(150))

sea_ice = sea_ice_simulation(grid; bottom_heat_boundary_condition, dynamics, advection=WENO(order=7))

set!(sea_ice.model.ice_thickness,     SI_meta_init; inpainting=nothing)
set!(sea_ice.model.ice_concentration, SC_meta_init; inpainting=nothing)

#####
##### Interface fluxes
#####

radiation = Radiation(sea_ice_albedo=0.7)
arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
arctic = Simulation(arctic, Δt=10, stop_time=30days)

ArcticOcean.arctic_outputs!(arctic, "EVP-rheology/")

# And add it as a callback to the simulation.
add_callback!(arctic, ArcticOcean.progress, IterationInterval(1))

# And add it as a callback to the simulation.
wizard = TimeStepWizard(cfl=0.7, max_change=1.001, max_Δt=10minutes)

function sea_ice_cell_advection_timescale(grid, velocities)
    u, v = velocities
    τ = KernelFunctionOperation{Center, Center, Center}(cell_advection_timescaleᶜᶜ, grid, u, v)
    return minimum(τ)
end

@inline _inverse_timescale(i, j, k, Δ⁻¹, U, topo) = @inbounds abs(U[i, j, k]) * Δ⁻¹
@inline _inverse_timescale(i, j, k, Δ⁻¹, U, topo::Flat) = 0

@inline function cell_advection_timescaleᶜᶜ(i, j, k, grid::AbstractGrid{FT, TX, TY}, u, v) where {FT, TX, TY}
    Δx⁻¹ = Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    Δy⁻¹ = Δy⁻¹ᶜᶠᶜ(i, j, k, grid)

    inverse_timescale_x = _inverse_timescale(i, j, k, Δx⁻¹, u, TX())
    inverse_timescale_y = _inverse_timescale(i, j, k, Δy⁻¹, v, TY())

    inverse_timescale = inverse_timescale_x + inverse_timescale_y

    return 1 / inverse_timescale
end

function add_wizard!(sim)
     wizard(sim.model.ocean)
     sea_ice = sim.model.sea_ice
     Δti = 0.1 * sea_ice_cell_advection_timescale(sea_ice.model.grid, sea_ice.model.velocities)
     @info "Wizard says: ocean Δt: $(ocean.Δt), sea ice Δt: $(Δti)"
     sim.Δt = min(ocean.Δt, Δti)
end

arctic.callbacks[:wizard] = Callback(add_wizard!, IterationInterval(1))

run!(arctic)
