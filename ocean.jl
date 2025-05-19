using ClimaOcean.DataWrangling: NearestNeighborInpainting

# A very diffusive ocean
momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

horizontal_Ri_filter = Oceananigans.TurbulenceClosures.FivePointHorizontalFilter()
free_surface = SplitExplicitFreeSurface(grid; substeps=100)
closure = Oceananigans.TurbulenceClosures.RiBasedVerticalDiffusivity(; horizontal_Ri_filter) # ClimaOcean.OceanSimulations.default_ocean_closure()

forcing = (T=RT, S=RS)

timestepper = :SplitRungeKutta3

ocean = ocean_simulation(grid;
                         momentum_advection,
                         tracer_advection,
                         timestepper,
                         free_surface,
                         forcing,
                         closure)

set!(ocean.model, T=T_meta_init, S=S_meta_init) 
