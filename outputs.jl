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

function arctic_outputs!(simulation)

    # Extract the sea ice and ocean models from the simulation
    sea_ice = simulation.model.sea_ice
    ocean   = simulation.model.ocean

    # Define the variables to be saved
    h = sea_ice.model.ice_thickness
    ℵ = sea_ice.model.ice_concentration
    u = sea_ice.model.velocities.u
    v = sea_ice.model.velocities.v
    
    # Fluxes
    Tu = ocean.model.interfaces.atmosphere_sea_ice_interface.temperature
    Qˡ = ocean.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
    Qˢ = ocean.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
    Qⁱ = ocean.model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat
    Qᶠ = ocean.model.interfaces.sea_ice_ocean_interface.fluxes.frazil_heat
    Qᵗ = ocean.model.interfaces.net_fluxes.sea_ice_top.heat
    Qᴮ = ocean.model.interfaces.net_fluxes.sea_ice_bottom.heat
    τx = ocean.model.interfaces.net_fluxes.sea_ice_top.u
    τy = ocean.model.interfaces.net_fluxes.sea_ice_top.v

    # Output writers for variables and averages
    simulation.output_writers[:vars] = JLD2Writer(sea_ice.model, (; h, ℵ, u, v, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ, τx, τy),
                                                  filename = "sea_ice_quantities.jld2",
                                                  schedule = IterationInterval(100),
                                                  overwrite_existing=true)

    simulation.output_writers[:averages] = JLD2Writer(sea_ice.model, (; h, ℵ, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ, u, v, τx, τy),
                                                      filename = "averaged_sea_ice_quantities.jld2",
                                                      schedule = AveragedTimeInterval(1days),
                                                      overwrite_existing=true)

    return nothing
end