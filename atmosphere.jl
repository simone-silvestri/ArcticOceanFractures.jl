using ClimaOcean.JRA55
using Oceananigans.OutputReaders: Cyclical

atmosphere = JRA55PrescribedAtmosphere(arch; time_indexing = Cyclical(), backend=JRA55NetCDFBackend(40), dir="./data") #, start_date, end_date, dir="./data") #, dataset=MultiYearJRA55(), start_date, end_date, dir="./data")
