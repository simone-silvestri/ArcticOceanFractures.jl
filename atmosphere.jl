using ClimaOcean.JRA55
using Oceananigans.OutputReaders: Cyclical

atmosphere = JRA55PrescribedAtmosphere(arch; dataset=MultiYearJRA55(), start_date, end_date, dir="./data")
