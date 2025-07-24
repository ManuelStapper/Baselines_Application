# Baseline Model Selection for Epidemiological Forecasting

This repository contains code and data processing scripts for the paper ["Mind the Baseline: The Hidden Impact of Reference Model Selection on Forecast Assessment"](https://github.com/epiforecasts/Baselines)

## Data Availability

All input data files, processed datasets, and analysis results are available on [Zenodo](https://doi.org/10.5281/zenodo.16407890)

## Data Sources

### Influenza Data
- [**FluSight Forecast Hub**](https://github.com/cdcepi/FluSight-forecast-hub): Forecast data for seasons 2023/24 and 2024/25
- [**FluSight Archive**](https://github.com/cdcepi/Flusight-forecast-data): Historical data (2021/22, 2022/23 seasons)
- **Target**: Weekly incident hospitalizations for US states
- **Truth Data**: Hospital admission counts from CDC surveillance

### COVID-19 Data
- [**COVID-19 Forecast Hub**](https://github.com/reichlab/covid19-forecast-hub): Forecast data
- **Target**: Weekly incident hospitalizations for US states
- **Truth Data**: Hospitalizations from publicly available sources via [covidHubUtils](https://github.com/reichlab/covidHubUtils)

## Requirements

### Julia Packages
Standard packages: CSV, DataFrames, Dates, Plots, StatsBase, Optim, Distributions, BoxCox

Visualization: Shapefile, GeoDataFrames, LibGEOS, ColorSchemes

### Custom Package
This analysis requires the BaselineModels.jl package:
```julia
using Pkg
Pkg.add(url="https://github.com/ManuelStapper/BaselineModels.jl", rev="afaafc8714239add0986243c44b58afad1765dbd")
```

## Code files
- **`DataPrep.jl`**: Cleans and filters raw forecast data from forecast hubs. Results are stored on [Zenodo](https://doi.org/10.5281/zenodo.16407890).
- **`InputPrep.jl`**: Converts cleaned data into structured forecast objects.
- **`CreateBaselineForecasts.jl`**: Runs forecasts for all baseline models and stores forecasts in iBase.csv and cBase.csv (see [Zenodo](https://doi.org/10.5281/zenodo.16407890))
- **`PlottingPrep.jl`**: Prepare shapefiles for plots
- **`Plotting.jl`**: Functions for plotting maps using above shapefiles
- **`Results.jl`**: Main analysis script




