# Set path here
# path = 
cd(path)
using CSV, DataFrames, Dates, Plots

# We need two main data sets
# Influenza & COVID

#################
### Influenza ###
#################

# There are two repositories for forecasts:
# FluSight 2023/2024 (?) 2024/2025
# https://github.com/cdcepi/FluSight-forecast-hub
# Is the currently used (also prev season?)
# FluSight 2021/2022 2022/2023
# https://github.com/cdcepi/Flusight-forecast-data

# Save only point and quantile forecasts
# Data format

# Read forecast model names from folders
cd(path*"/Influenza/FluSightNew/model-output")
models = readdir()
models = models[(m -> length(split(m, ".")) == 1).(models)]

# Each file has columns:
# - reference_date:     When are forecasts conducted?
# - target              Should be "wk inc flu hosp"
# - horizon             Values between -1 and 3 (what do -1 and 0 mean?)
# - target_end_date     Week to be forecast
# - location            Numeric or "US" referencing states / US
# - output_type         Quantile, point, pmf... (What is the forecast?) Only quantile and point 
# - output_type_id      If quantile, gives the percentage 
# - value               Forecast (quantile or point)

# Restrict locations to 1, ..., 56 or "US"
# "US" -> 0, exclude terretories

function stringToInteger(x::String)
    if x == "US"
        return 0
    end
    parse(Int64, x)
end

function cleanData(x::DataFrame, modelName::String)
    indType = (x.output_type .== "quantile")
    indTarget = x.target .== "wk inc flu hosp"
    if typeof(x.location) == Vector{Int64}
        stateInt = x.location
    else
        stateInt = stringToInteger.(Vector{String}(x.location))
    end
    
    indLocation = stateInt .<= 56
    indKeep = indType .& indTarget .& indLocation
    x = x[indKeep, :]
    x.location = stateInt[indKeep]
    x.model = fill(modelName, nrow(x))
    return x
end


# For each model folder browse through all csv files
# Problem: There is one file with too many rows causing an error
# Solution: Removed forecasts that are not needed (not influenza hospitalisation) "by hand"

dfVector = Vector{DataFrame}(undef, 0)

for m = 1:length(models)
    cd(path* "/Influenza/FluSightNew/model-output/"*models[m])
    files = readdir()
    files = files[(f -> split(f, ".")[2] == "csv").(files)]
    for i = 1:length(files)
        temp = CSV.read(files[i], DataFrame)
        temp = cleanData(temp, models[m])
        push!(dfVector, temp)
    end
end

influenzaNew = vcat(dfVector...)


# "Old" data

# Each file has columns:
# - forecast_date:     When are forecasts conducted?
# - target              Should be "n wk inc flu hosp" for n = 1, 2, 3, 4
# - target_end_date     Week to be forecast
# - location            Numeric or "US" referencing states / US
# - type                Quantile, point, pmf... (What is the forecast?) Only quantile and point 
# - quantile            If quantile, gives the percentage, otherwise "NA"
# - value               Forecast (quantile or point)

# Translate to above notation of new data


function cleanDataOld(x::DataFrame, modelName::String)
    rename!(x, :forecast_date => :reference_date, :type => :output_type, :quantile => :output_type_id)
    indType = (x.output_type .== "quantile")
    indTarget = (tt -> tt in string.(1:4) .* " wk ahead inc flu hosp").(x.target)
    x = x[indType .& indTarget, :]

    if typeof(x.location) == Vector{Int64}
        stateInt = x.location
    else
        stateInt = stringToInteger.(Vector{String}(x.location))
    end
    
    indLocation = stateInt .<= 56
    x = x[indLocation, :]
    x.location = stateInt[indLocation]
    x.horizon = (tt -> parse(Int64, split(tt, " ")[1])).(x.target)
    x.target = fill("wk inc flu hosp", nrow(x))
    x.model = fill(modelName, nrow(x))
    
    x = x[:, [1, 2, 8, 3, 4, 5, 6, 7, 9]]
    return x
end

cd(path*"/Influenza/FluSightOld/data-forecasts")
models = readdir()
models = models[(m -> length(split(m, ".")) == 1).(models)]

dfVector = Vector{DataFrame}(undef, 0)

for m = 1:length(models)
    cd(path* "/Influenza/FluSightOld/data-forecasts/"*models[m])
    files = readdir()
    files = files[(f -> split(f, ".")[2] == "csv").(files)]
    for i = 1:length(files)
        temp = CSV.read(files[i], DataFrame)
        temp = cleanDataOld(temp, models[m])
        push!(dfVector, temp)
    end
end

function replaceDate(x)
    if typeof(x) == Date
        return x
    end
    for splitChar = ["-", "/", "."]
        temp = split(x, splitChar)
        if length(temp) == 3
            if parse(Int64, temp[1]) < 1000
                temp = reverse(temp)
            end
            return Date(parse(Int64, temp[1]), parse(Int64, temp[2]), parse(Int64, temp[3]))
        end
    end
end

influenzaOld = vcat(dfVector...)

influenzaOld = influenzaOld[.!ismissing.(influenzaOld.value), :]
influenzaOld.value = Vector{Float64}(influenzaOld.value)

indString = typeof.(influenzaNew.reference_date) .== String
influenzaNew.reference_date[indString] .= Date(2024, 11, 30)
influenzaNew.reference_date = Vector{Date}(influenzaNew.reference_date)

influenzaNew.target_end_date = replaceDate.(influenzaNew.target_end_date)

function replaceOutputTypeId(x)
    if typeof(x) == Float64
        return x
    end
    return parse(Float64, x)
end

influenzaNew.output_type_id = replaceOutputTypeId.(influenzaNew.output_type_id)
influenzaOld.output_type_id = replaceOutputTypeId.(influenzaOld.output_type_id)

influenzaNew.horizon[typeof.(influenzaNew.horizon) .== String] = parse.(Int64, influenzaNew.horizon[typeof.(influenzaNew.horizon) .== String])
influenzaNew.horizon = Vector{Int64}(influenzaNew.horizon)


# Done with forecast data
influenza = [influenzaNew; influenzaOld]

influenza.model = uppercase.(influenza.model)
# No duplicate rows when only considering target_end_date, horizon, location, output_type_id and model

# Remove constant variables
influenza = influenza[:, [1, 3, 4, 5, 7, 8, 9]]

models = sort(unique(influenza.model))

modelID = (x -> findfirst(x .== models)).(influenza.model)
influenza.model = modelID

influenza.reference_date = influenza.target_end_date .- influenza.horizon .* Day(7)
sort!(influenza, ["reference_date", "location", "model", "horizon", "output_type_id"])

forecastID = ones(Int64, nrow(influenza))
counter = 1
for i = 2:nrow(influenza)
    if influenza[i-1, [1, 4, 7]] != influenza[i, [1, 4, 7]]
        counter += 1
    end
    forecastID[i] = counter
end
influenza.forecastID = forecastID

cd(path*"/Influenza")
# CSV.write("Influenza.csv", influenza)
# CSV.write("InfluenzaModelID.csv", DataFrame(ID = 1:length(models), model = models))


# Influenza truth data

cd(path * "/Influenza/FluSightNew/target-data")
truthNew = CSV.read("target-hospital-admissions.csv", DataFrame)

truthNew.location = stringToInteger.(truthNew.location)
truthNew = truthNew[truthNew.location .<= 56, :]

truthNew = truthNew[:, [1, 2, 4]]

truthNew.value = Vector{Union{String, Missing, Int64}}(truthNew.value)
truthNew.value[truthNew.value .== "NA"] .= missing
truthNew.value[.!ismissing.(truthNew.value)] = parse.(Int64, truthNew.value[.!ismissing.(truthNew.value)])
truthNew.value = Vector{Union{Int64, Missing}}(truthNew.value)

cd(path * "/Influenza/FluSightOld/data-truth")
truthOld = CSV.read("truth-Incident Hospitalizations.csv", DataFrame)

truthOld.location = stringToInteger.(truthOld.location)
truthOld = truthOld[truthOld.location .<= 56, :]

truthOld = truthOld[:, [1, 2, 4]]

truth = [truthNew; truthOld]

cd(path*"/Influenza")
# CSV.write("InfluenzaTruth.csv", truth)


#############
### COVID ###
#############

# COVID forecast data
# COVID truth data

# Massive data set
# Filter the forecast data and delete folders
# Afterwards add truth data

# Function to filter out data that has any useful forecast
function cleanCOVIDinitial(x::DataFrame, model::String)
    keepTarget = (z -> z[3:end] .== "wk ahead inc case").(x.target)
    keepQuantile = x.type .== "quantile"

    x = x[keepTarget .& keepQuantile, :]
    if sum(keepTarget .& keepQuantile) == 0
        return DataFrame(forecast_date = Date[], target_end_date = Date[], horizon = Int64[], location = Int64[], quantile = Float64[], value = Float64[], model = String[])
    else
        x.model = fill(model, nrow(x))
        return x
    end
end

# Read in forecast data first to see which truth data is needed
cd(path * "/COVID/covid19-forecast-hub/data-processed")
models = readdir()
models = models[(m -> length(split(m, ".")) == 1).(models)]

# For each model, aggregate forecast data and delete folder

# Notes:
# Problem with model "JHUAPL-SLPHospEns": daily forecasts makes forecast files too large to be read in
# --> Remove model
# Maybe there are other models with the same problem

for m = 1:length(models)
    dfVector = Vector{DataFrame}(undef, 0)
    cd(path * "/COVID/covid19-forecast-hub/data-processed/" * models[m])
    files = readdir()
    files = files[(f -> length(split(f, ".")) == 2).(files)]
    files = files[(f -> split(f, ".")[2] == "csv").(files)]
    for f = files
        temp = CSV.read(f, DataFrame)
        temp = cleanCOVIDinitial(temp, models[m])
        if nrow(temp) > 0
            push!(dfVector, temp)
        end
    end
    
    # If there is no useful forecast data, remove folder and move on
    if length(dfVector) == 0
        rm(path * "/COVID/covid19-forecast-hub/data-processed/" * models[m], recursive=true)
    else
        # If there are "valid" forecasts, create new folder with aggregated forecasts
        # and remove the raw data folder to save space
        dfVector = vcat(dfVector...)
        mkdir(path*"/COVID/ForecastsAggregated/"*models[m])
        CSV.write(path*"/COVID/ForecastsAggregated/"*models[m]*"/forecasts.csv", dfVector)
        rm(path * "/COVID/covid19-forecast-hub/data-processed/" * models[m], recursive=true)
    end    
    println(m)
end

# Filter out counties
cd(path*"/COVID/ForecastsAggregated")
models = readdir()
models = models[(m -> length(split(m, ".")) == 1).(models)]

mSkip = Int64[]

for m = 1:length(models)
    temp = CSV.read(path*"/COVID/ForecastsAggregated/"*models[m]*"/forecasts.csv", DataFrame)
    if typeof(temp.location) != Vector{Int64}
        temp.location = stringToInteger.(temp.location)
        temp = temp[temp.location .<= 56, :]
        CSV.write(path*"/COVID/ForecastsAggregated/"*models[m]*"/forecasts.csv", temp)
        println(m)
    else
        push!(mSkip, m)
    end
end

for m = mSkip
    rm(path*"/COVID/ForecastsAggregated/"*models[m], recursive=true)
end

# Now do the cleaning and aggregate all forecasts
# Goal: Have one file with all relevant forecasts that can later be structured in "forecast" object types

cd(path*"/COVID/ForecastsAggregated")
models = readdir()
models = models[(m -> length(split(m, ".")) == 1).(models)]
# 58 models remaining

dfVector = Vector{DataFrame}(undef, 0)

# Target format:
# forecast_date, location, horizon, model, quantiles

for m = 1:length(models)
    temp = CSV.read(path*"/COVID/ForecastsAggregated/"*models[m]*"/forecasts.csv", DataFrame)
    temp = temp[temp.value .!= "NULL", :]
    if (typeof(temp.value) != Vector{Int64}) & (typeof(temp.value) != Vector{Float64})
        temp.value = parse.(Float64, temp.value)
    else
        temp.value = Vector{Float64}(temp.value)
    end
    temp.horizon = parse.(Int64, (x -> split(x, " ")[1]).(temp.target))
    temp.forecast_date = temp.target_end_date .- temp.horizon .* Day(7)
    out = DataFrame(forecast_date = temp.forecast_date,
                    horizon = temp.horizon,
                    location = temp.location,
                    quantile = temp.quantile,
                    value = temp.value,
                    model = temp.model)
    push!(dfVector, out)
end

covid = vcat(dfVector...)

# Remove observations with horizon > 4

function findQuantiles(df::DataFrame, q::Vector{Float64})
    out = zeros(length(q))
    for i = 1:length(q)
        keep = df.quantile .== q[i]
        if sum(keep) == 0
            out[i] = -1
        else
            out[i] = df.value[findfirst(keep)]
        end
    end
    out
end

covid = covid[covid.horizon .<= 4, :]
models = sort(unique(covid.model))
quantiles = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]

dfVector = Vector{DataFrame}(undef, 0)

for m = 1:length(models)
    temp = covid[covid.model .== models[m], :]
    # Unique combionations to be added:
    tempUni = unique(temp[:, ["forecast_date", "location", "horizon"]])
    vals = zeros(nrow(tempUni), 7)
    for i = 1:nrow(tempUni)
        ind1 = (temp.forecast_date .== tempUni.forecast_date[i])
        ind2 = (temp.location .== tempUni.location[i])
        ind3 = (temp.horizon .== tempUni.horizon[i])
        temp2 = temp[ind1 .& ind2 .& ind3, :]
        vals[i, :] = findQuantiles(temp2, quantiles)
    end

    vals_df = DataFrame(vals, :auto)
    rename!(vals_df, ["q025", "q100", "q250", "q500", "q750", "q900", "q975"])
    tempUni = [tempUni vals_df]
    tempUni.model = fill(models[m], nrow(tempUni))
    push!(dfVector, tempUni)
    println(m)
end

covid2 = vcat(dfVector...)
# Remove forecasts with missings

indKeep = all(Matrix{Float64}(covid2[:, 4:10]) .> -1, dims = 2)[:, 1]
covid2 = covid2[indKeep, :]

CSV.write(path*"/COVID/covid.csv", covid2)

# For truth data covid, download and install git large files
# truth data fetched from there and saved in a separate folder

truth = CSV.read(path*"/COVID/TruthGitLFS/truth-Incident Cases.csv", DataFrame)
truth

truth.location = stringToInteger.(truth.location)
truth = truth[truth.location .<= 56, :]
# Aggregate weeks (sun - sat) with saturday being the "target_date"

truth = sort!(truth, ["location", "date"])

startSeq = collect(Date(2020, 01, 05):Day(7):Date(2024, 04, 13))
endSeq = startSeq .+ Day(6)

edVec = Date[]
lVec = Int64[]
valVec = Int64[]

for l = sort(unique(truth.location))
    for i = 1:length(startSeq)
        # If there are less than 7 dates available, treat as missing
        ind = (truth.location .== l) .& (startSeq[i] .<= truth.date .<= endSeq[i])
        if sum(ind) == 7
            push!(valVec, sum(truth.value[ind]))
            push!(lVec, l)
            push!(edVec, endSeq[i])
        end
    end
end


truth = DataFrame(location = lVec, end_date = edVec, value = valVec)

CSV.write(path*"/COVID/truth.csv", truth)




############
### Misc ###
############

# Using the data, we can create overviews of availability for each model

models = sort(unique(influenza.model))
target_dates = sort(unique(influenza.target_end_date))
ref_dates = sort(unique(influenza.reference_date))
locations = sort(unique(influenza.location))
q = sort(unique(influenza.output_type_id))

# 86 models, 110 dates to be forecast, 52 regions, 4 horizons being forecast (which differs between new and old!)
# and also 23 different quantiles

unique(influenza.reference_date)
nObs = (m -> nrow(influenza[influenza.model .== m, :])).(models)

histogram(nObs / (110 * 52 * 4 * 23)*100, nbins = 20)

avail = zeros(Int64, length(models), 5)
cutoffDates = [Date(2022, 07, 31), Date(2023, 07, 31), Date(2024, 07, 31)]

# Number of dates to be forecast by season:
# 27, 36, 35, 12

for m = 1:length(models)
    temp = influenza[influenza.model .== models[m], :]
    avail[m, 1] = length(unique(temp.location))
    tempDates = sort(unique(temp.target_end_date))
    tempSeason = (tt -> sum(tt .> cutoffDates) + 1).(tempDates)
    avail[m, 2] = sum(tempSeason .== 1)
    avail[m, 3] = sum(tempSeason .== 2)
    avail[m, 4] = sum(tempSeason .== 3)
    avail[m, 5] = sum(tempSeason .== 4)
end