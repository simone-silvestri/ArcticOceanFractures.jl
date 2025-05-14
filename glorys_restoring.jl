using PythonCall
using ClimaOcean
using ClimaOcean.Copernicus
using ClimaOcean.ECCO
using ClimaOcean.DataWrangling
using Oceananigans
using Oceananigans.OutputReaders

bounding_box = DataWrangling.BoundingBox(longitude=(-180.0, 180.0), latitude=(0.0, 90.0), z=(-6000, 0))

T_meta_init = Metadatum(:temperature; dataset=GLORYSDaily(), date=start_date, bounding_box, dir="data/") 
S_meta_init = Metadatum(:salinity;    dataset=GLORYSDaily(), date=start_date, bounding_box, dir="data/") 

SI_meta_init = Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly(), date=start_date, dir="data/")
SC_meta_init = Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly(), date=start_date, dir="data/")

# Restore quite strongly between 40 and 42.5 north then decrease linearly between 42.5 and 45 north
function mask(λ, φ, z, t) 
    m = ifelse(φ < 85//2, one(φ),
        ifelse(φ < 45,    1 - (φ - 85//2) / (45 - 85//2), 
                          zero(φ)))
    return m
end

T_meta_rest = Metadata(:temperature; dataset=ECCO4Monthly(), start_date, end_date, dir="data/") 
S_meta_rest = Metadata(:salinity;    dataset=ECCO4Monthly(), start_date, end_date, dir="data/") 

kwargs = (; rate = 1/10days, mask, time_indexing=Linear(), time_indices_in_memory=5)

RT = ClimaOcean.DatasetRestoring(T_meta_rest, arch; kwargs...)
RS = ClimaOcean.DatasetRestoring(S_meta_rest, arch; kwargs...)
