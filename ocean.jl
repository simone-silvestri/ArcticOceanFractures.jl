
# A very diffusive ocean
momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; substeps=100)
closure = ClimaOcean.OceanSimulations.default_ocean_closure()

forcing = (T=RT, S=RS)

ocean = ocean_simulation(grid;
                         momentum_advection,
                         tracer_advection,
                         free_surface,
                         forcing,
                         closure)

set!(ocean.model; T=first(T_meta_init), S=first(S_meta_init))