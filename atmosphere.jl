using ClimaOcean.JRA55

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(40), dataset=MultiYearJRA55(), start_date, end_date, dir="./data")
