#########################################################
### Transform cleaned data into managable data frames ###
#########################################################

# Set path here
# path = 
cd(path)
using CSV, DataFrames, Dates, Plots, BaselineModels

#################
### Influenza ###
#################

cd(path*"/Influenza")
influenza = CSV.read("Influenza.csv", DataFrame)
influenzaModelID = CSV.read("InfluenzaModelID.csv", DataFrame)
influenzaTruth = CSV.read("InfluenzaTruth.csv", DataFrame)

# Create forecast object for each
# Reference date, model, location

# Create a cleaned up data frame with reference_date x model x location and store forecast

q = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99]

function dfToForecast(x::DataFrame)
    ### Takes a data frame and translates to one forecast object

    horizon = sort(unique(x.horizon))

    # initialise intervals and median vector
    intervals = forecastInterval[]
    med = Float64[]

    for h = horizon
        # Find part of data frame for horizon h
        tempH = x[x.horizon .== h, :]
        push!(med, tempH.value[tempH.output_type_id .== 0.5][1])
        # Filter out asymmetric intervals
        quantiles = tempH.output_type_id
        q1 = quantiles[quantiles .< 0.5]
        q2 = round.(1 .- quantiles[quantiles .> 0.5], digits = 3)
        qKeep = intersect(q1, q2)
        qKeep = round.(sort(unique([qKeep; 1 .- qKeep])), digits = 3)
        qKeep = (qq -> qq in qKeep).(quantiles)
        tempH = tempH[qKeep, :]

        quantiles = tempH.output_type_id
        # Translate quantiles to α level
        α = 2 .* quantiles[quantiles .< 0.5]
        # Read lower and upper bounds of intervals
        L = tempH.value[quantiles .< 0.5]
        U = reverse(tempH.value[quantiles .> 0.5])

        if length(L) != length(U)
            println(tempH)
        end
        push!(intervals, forecastInterval(α, L, U))
        # Get median under assumption that there is a median prediction
        # push!(med, tempH.value[sum(quantiles .< 0.5) + 1])
    end
    forecast(horizon, median = med, interval = intervals)
end

# Function to translate the whole data frame to forecast objects
# It has a variable forecastID, that separates the forecast objects foe speed up
# It returns a data frame with forecast objects (and reference date, model, locations)
function prep1(influenza)
    # Initialise variables of output data frame
    forecastVec = forecast[] # For forecast
    rdVec = Date[] # For reference date
    mVec = Int64[] # For model
    lVec = Int64[] # For location
    
    # Index variables for relevant rows of influenza
    iL = 1
    iU = 1
    for i = 2:nrow(influenza)-1
        if influenza.forecastID[i] == influenza.forecastID[i-1]
            # Not yet the end of forecast block
            iU += 1
        else
            # If the end of the forecast block is found 
            # Here i is the start of a new block
            temp = influenza[iL:iU, :]
            # Reference date, model and location the same for all rows of the block
            push!(rdVec, temp.reference_date[1])
            push!(mVec, temp.model[1])
            push!(lVec, temp.location[1])
            # Translate the block to forecast object
            push!(forecastVec, dfToForecast(temp))
            iL = i
            iU = i
        end
    end
    # Final block:    
    temp = influenza[iL:end, :]
    push!(rdVec, temp.reference_date[1])
    push!(mVec, temp.model[1])
    push!(lVec, temp.location[1])
    push!(forecastVec, dfToForecast(temp))

    return rdVec, mVec, lVec, forecastVec
end

rdVec, mVec, lVec, forecastVec = prep1(influenza)

function checkMissing(x::forecast)
    sum((i -> sum((i.l .== -1) .| (i.u .== -1))).(x.interval) .> 0)
end
# None missing
# nMissing = checkMissing.(forecastVec)

# Fill forecast objects with truth data
# Resolve missings and duplicates in influenzaTruth
dropmissing!(influenzaTruth, :value)
influenzaTruth = combine(groupby(influenzaTruth, [:date, :location])) do sdf
    sdf[findmax(sdf.value)[2], :]
end

### Fill in truth data, iterate over all forecasts
### (Quite time consuming)
for i = 1:length(forecastVec)
    # Find out target dates of forecasts
    dateSeq = rdVec[i] .+ forecastVec[i].horizon .* Day(7)

    # Filter the corresponding rows in truth data
    # ind = findall((influenzaTruth.location .== lVec[i]) .& (dd -> dd in dateSeq).(influenzaTruth.date))
    ind = findall((influenzaTruth.location .== lVec[i]) .& (influenzaTruth.date .<= dateSeq[end]) .& (influenzaTruth.date .>= dateSeq[1]))
    truthTemp = fill(-1, length(forecastVec[i].horizon))

    for j in 1:length(ind)
        truthTemp[dateSeq .== influenzaTruth.date[ind[j]]] .= influenzaTruth.value[ind][j]
    end
    forecastVec[i].truth = truthTemp
end

ifc = DataFrame(date = rdVec, model = mVec, location = lVec, forecast = forecastVec)

#############
### COVID ###
#############

cd(path*"/COVID")
covid = CSV.read("covid.csv", DataFrame)
covidTruth = CSV.read("truth.csv", DataFrame)

covidUni = unique(covid[:, ["forecast_date", "location", "model"]])

function prep2(covidUni, covid)
    rdVec = Date[]
    mVec = String[]
    lVec = Int64[]
    forecastVec = forecast[]

    currModel = ""
    covidTemp = covid[covid.model .== currModel, :]

    for i = 1:nrow(covidUni)
        if covidUni.model[i] != currModel
            currModel = covidUni.model[i]
            covidTemp = covid[covid.model .== currModel, :]
        end
        ind1 = covidTemp.location .== covidUni.location[i]
        ind2 = covidTemp.forecast_date .== covidUni.forecast_date[i]
        ind3 = covidTemp.model .== covidUni.model[i]
        temp = covidTemp[ind1 .& ind2 .& ind3, :]

        temp = temp[sortperm(temp.horizon), :]
        
        intervals = (i -> forecastInterval([0.05, 0.2, 0.5], Vector{Float64}(temp[i, 4:6]), Vector{Float64}(temp[i, 10:-1:8]))).(1:nrow(temp))
        med = temp.q500
        fc = forecast(temp.horizon, median = med, interval = intervals)

        push!(rdVec, covidUni.forecast_date[i])
        push!(mVec, covidUni.model[i])
        push!(lVec, covidUni.location[i])
        push!(forecastVec, fc)
    end
    return rdVec, mVec, lVec, forecastVec
end

rdVec, mVec, lVec, forecastVec = prep2(covidUni, covid)


# Add truth data to forecasts

function getTruth(date, location)
    ind = (covidTruth.location .== location) .& (covidTruth.end_date .== date)
    if sum(ind) != 1
        return -1
    else
        return covidTruth.value[ind][1]
    end
end

for i = 1:length(rdVec)
    horizons = forecastVec[i].horizon
    targets = rdVec[i] .+ Day(7) .* horizons
    truth = (tt -> getTruth(tt, lVec[i])).(targets)
    forecastVec[i].truth = truth
end

cfc = DataFrame(date = rdVec, model = mVec, location = lVec, forecast = forecastVec)

models = sort(unique(covid.model))
covidModelID = DataFrame(ID = collect(1:length(models)), model = models)
# CSV.write("covidModelID.csv", covidModelID)
cfc.model = (mm -> findfirst(mm .== models)).(mVec)

sort!(ifc, ["model", "location", "date"])
sort!(cfc, ["model", "location", "date"])

sort!(influenzaTruth, ["location", "date"])
rename!(covidTruth, ["location", "date", "value"])
sort!(covidTruth, ["location", "date"])
