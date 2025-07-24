#############################################
### Creating forecasts of baseline models ###
#############################################

using CountTimeSeries


using CSV, DataFrames, Dates, Plots
using BaselineModels
# using StatsBase, Random, Optim

# Set path here
# path = 
include(path*"/InputPrep.jl")
cd(path)


# Forecast data needs to be cleaned further:
# Only include positive horizon forecasts and forecasts where truth data is available
keep = fill(true, nrow(ifc))

for i = 1:nrow(ifc)
    fc = ifc.forecast[i]
    ind = findall(fc.horizon .> 0)
    ind = setdiff(ind, findall(fc.truth .== -1))
    if length(ind) > 0
        ifc.forecast[i] = forecast(fc.horizon[ind], interval = fc.interval[ind], truth = fc.truth[ind], median = fc.median[ind])
    else
        keep[i] = false
    end
end
ifc = ifc[keep, :]

dates = sort(unique(ifc.date)) # Reference dates
locations = sort(unique(ifc.location))

# Apply good models here, will later be contrasted with
# bad baseline models that look like a good choice?
models = Vector{Baseline}(undef, 19)
models[1] = constantModel(true)
models[2] = marginalModel(52, true)
models[3] = kdeModel(true)
models[4] = lsdModel(52, 5, true)
models[5] = olsModel(5, 1, true)
models[6] = idsModel(3, true)
models[7] = armaModel(1, 1, m = 52, trend = false, isPos = true)
models[8] = inarchModel(1, 52, 8, true)
models[9] = etsModel(error = "A", trend = "N", season = "N")
models[10] = stlModel(p = 1, o = 1)

models[11] = constantModel(false)
models[12] = marginalModel(52, false)
models[13] = kdeModel(false)
models[14] = lsdModel(52, 5, false)
models[15] = olsModel(5, 1, false)
models[16] = idsModel(3, false)
models[17] = armaModel(1, 1, m = 52, trend = false, isPos = false)
models[18] = etsModel(error = "A", trend = "N", season = "N")
models[19] = stlModel(p = 1, o = 1)

quantiles = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99]
nD = length(dates)
nM = length(models)
nL = length(locations)


iBase = DataFrame(date = repeat(repeat(dates, outer = nL), nM),
                  model = repeat(1:nM, inner = nD*nL),
                  location = repeat(repeat(locations, inner = nD), nM),
                  forecast = Vector{forecast}(undef, nM*nD*nL))

# Filter out the respective truth data
temp = influenzaTruth[influenzaTruth.location .== iBase.location[1], :]
x = temp.value
xl = log.(x .+ 1)
tEnd = findlast(temp.date .<= iBase.date[1])

h = minimum([4, length(x) - tEnd])
res = fit(x[1:tEnd], models[iBase.model[1]])
fc = BaselineModels.predict(res, h, quantiles)
fc.truth = x[tEnd+1:tEnd+h]
iBase.forecast[1] = fc

for i = 2:nrow(iBase)
    if iBase.model[i] != iBase.model[i-1]
        println(now())
        println("Finished model " * string(iBase.model[i-1]))
        println("")
    end
    if iBase.location[i] != iBase.location[i-1]
        temp = influenzaTruth[influenzaTruth.location .== iBase.location[i], :]
        x = temp.value
        xl = log.(x .+ 1)
    end
    tEnd = findlast(temp.date .<= iBase.date[i])
    h = minimum([4, length(x) - tEnd])

    if (iBase.model[i] == 9) | (iBase.model[i] == 10)
        xHat, μ, est = preFilter(x[1:tEnd], 52, 2)
        res = fit(xHat, models[iBase.model[i]])
        fc = predict(res, h, quantiles)
        fc2 = postFilter(x[1:tEnd], fc, 52, est)
        fc2.truth = x[tEnd+1:tEnd+h]
        iBase.forecast[i] = fc2
    elseif (iBase.model[i] == 18) | (iBase.model == 19)
        xHat, μ, est = preFilter(xl[1:tEnd], 52, 2)
        res = fit(xHat, models[iBase.model[i]])
        fc = predict(res, h, quantiles)
        fc2 = postFilter(xl[1:tEnd], fc, 52, est)
        fc2.mean = Float64[]
        fc2.median = exp.(fc2.median) .- 1
        fc2.interval = (xx -> forecastInterval(xx.α, exp.(xx.l), exp.(xx.u))).(fc2.interval)
        fc2.truth = x[tEnd+1:tEnd+h]
        iBase.forecast[i] = fc2
    elseif iBase.model[i] <= 9
        res = fit(x[1:tEnd], models[iBase.model[i]])
        fc = BaselineModels.predict(res, h, quantiles)
        fc.truth = x[tEnd+1:tEnd+h]
        iBase.forecast[i] = fc
    else
        res = fit(xl[1:tEnd], models[iBase.model[i]])
        fc = BaselineModels.predict(res, h, quantiles)
        fc.mean = Float64[]
        fc.median = exp.(fc.median) .- 1
        fc.interval = (xx -> forecastInterval(xx.α, exp.(xx.l) .- 1, exp.(xx.u) .- 1)).(fc.interval)
        fc.truth = x[tEnd+1:tEnd+h]
        iBase.forecast[i] = fc
    end
end

# COVID

keep = fill(true, nrow(cfc))

for i = 1:nrow(cfc)
    fc = cfc.forecast[i]
    ind = findall(fc.horizon .> 0)
    ind = setdiff(ind, findall(fc.truth .== -1))
    if length(ind) > 0
        cfc.forecast[i] = forecast(fc.horizon[ind], interval = fc.interval[ind], truth = fc.truth[ind], median = fc.median[ind])
    else
        keep[i] = false
    end
end
cfc = cfc[keep, :]

# Visiaul check: There are negative values in the truth data
# Simply flipping the sign should be okay

covidTruth.value = abs.(covidTruth.value)

# Now let's look at the time seriens
locations = sort(unique(cfc.location))
dates = sort(unique(cfc.date))

# Weird thing going on at location locations[11]
l = 11
d1 = covidTruth.date[covidTruth.location .== locations[l]]
t1 = covidTruth.value[covidTruth.location .== locations[l]]

covidTruth.value[(covidTruth.date .>= d1[112]) .& (covidTruth.location .== locations[11])] .= floor.(Int64, (t1[111:end-1] .+ t1[112:end])./2)

modelsC = Vector{Baseline}(undef, 19)
modelsC[1] = constantModel(true)
modelsC[2] = marginalModel(1000, true)
# LSD model without seasonality is simply a marginal model
modelsC[3] = marginalModel(5, true)
modelsC[4] = kdeModel(true)
modelsC[5] = olsModel(5, 1, true)
modelsC[6] = idsModel(3, true)
modelsC[7] = armaModel(1, 1, m = 0, trend = true, isPos = true)
modelsC[8] = inarchModel(1, 1, 0, true)
modelsC[9] = etsModel(error = "M", trend = "Ad", season = "N")
modelsC[10] = stlModel(p = 1, o = 1)

modelsC[11] = constantModel(false)
modelsC[12] = marginalModel(5, false)
modelsC[13] = marginalModel(1000, true)
modelsC[14] = kdeModel(false)
modelsC[15] = olsModel(5, 1, false)
modelsC[16] = idsModel(3, false)
modelsC[17] = armaModel(1, 1, m = 0, trend = true, isPos = false)
modelsC[18] = etsModel(error = "A", trend = "Ad", season = "N")
modelsC[19] = stlModel(p = 1, o = 1)

quantiles = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99]

nD = length(dates)
nM = length(modelsC)
nL = length(locations)

cBase = DataFrame(date = repeat(repeat(dates, outer = nL), nM),
                  model = repeat(1:nM, inner = nD*nL),
                  location = repeat(repeat(locations, inner = nD), nM),
                  forecast = Vector{forecast}(undef, nM*nD*nL))

# Filter out the respective truth data
temp = covidTruth[covidTruth.location .== cBase.location[1], :]
x = temp.value
xl = log.(x .+ 1)
tEnd = findlast(temp.date .<= cBase.date[1])

h = minimum([4, length(x) - tEnd])
res = fit(x[1:tEnd], modelsC[cBase.model[1]])
fc = predict(res, h, quantiles)
fc.truth = x[tEnd+1:tEnd+h]
cBase.forecast[1] = fc

for i = 2:nrow(cBase)
    if cBase.model[i] != cBase.model[i-1]
        println(now())
        println("Finished model " * string(cBase.model[i-1]))
        println("")
    end
    if cBase.location[i] != cBase.location[i-1]
        temp = covidTruth[covidTruth.location .== cBase.location[i], :]
        x = temp.value
        xl = log.(x .+ 1)
    end
    tEnd = findlast(temp.date .<= cBase.date[i])
    h = minimum([4, length(x) - tEnd])

    if cBase.model[i] <= 9
        res = fit(x[1:tEnd], modelsC[cBase.model[i]])
        fc = predict(res, h, quantiles)
        fc.truth = x[tEnd+1:tEnd+h]
        cBase.forecast[i] = fc
    else
        res = fit(xl[1:tEnd], modelsC[cBase.model[i]])
        fc = predict(res, h, quantiles)
        fc.mean = Float64[]
        fc.median = exp.(fc.median) .- 1
        fc.interval = (xx -> forecastInterval(xx.α, exp.(xx.l) .- 1, exp.(xx.u) .- 1)).(fc.interval)
        fc.truth = x[tEnd+1:tEnd+h]
        cBase.forecast[i] = fc
    end
end

# Double check if there are negatives in intervals
function makePos(x::Vector{T}) where {T <: Real}
    if minimum(x) < 0
        x[x .< 0] .= 0
    end
    x
end

for i = 1:nrow(iBase)
    for hh = 1:length(iBase.forecast[i].interval)
        iBase.forecast[i].interval[hh].l = makePos(iBase.forecast[i].interval[hh].l)
        iBase.forecast[i].interval[hh].u = makePos(iBase.forecast[i].interval[hh].u)
    end
    iBase.forecast[i].median = makePos(iBase.forecast[i].median)
end


for i = 1:nrow(cBase)
    for hh = 1:length(cBase.forecast[i].interval)
        cBase.forecast[i].interval[hh].l = makePos(cBase.forecast[i].interval[hh].l)
        cBase.forecast[i].interval[hh].u = makePos(cBase.forecast[i].interval[hh].u)
    end
    cBase.forecast[i].median = makePos(cBase.forecast[i].median)
end


# Saving the forecasts

function forecast2vec(x::forecast, h::Int64)
    [x.median[h]; x.truth[h]; x.interval[h].l; x.interval[h].u]
end

N = sum((x -> length(x.horizon)).(iBase.forecast))

iBaseMat = zeros(N, 24)
iBaseDat = fill(Date(2000, 1, 1), N)
iBaseMod = zeros(Int64, N)
iBaseLoc = zeros(Int64, N)
iBaseHor = zeros(Int64, N)

counter = 1
for i = 1:nrow(iBase)
    for h = iBase.forecast[i].horizon
        iBaseMat[counter, :] = forecast2vec(iBase.forecast[i], h)
        iBaseDat[counter] = iBase.date[i]
        iBaseMod[counter] = iBase.model[i]
        iBaseLoc[counter] = iBase.location[i]
        iBaseHor[counter] = h
        counter += 1
    end
end


iBaseDF = [DataFrame(date = iBaseDat, model = iBaseMod, location = iBaseLoc, horizon = iBaseHor) DataFrame(iBaseMat, :auto)]
# CSV.write("iBase.csv", iBaseDF)


N = sum((x -> length(x.horizon)).(cBase.forecast))

cBaseMat = zeros(N, 24)
cBaseDat = fill(Date(2000, 1, 1), N)
cBaseMod = zeros(Int64, N)
cBaseLoc = zeros(Int64, N)
cBaseHor = zeros(Int64, N)

counter = 1
for i = 1:nrow(cBase)
    for h = cBase.forecast[i].horizon
        cBaseMat[counter, :] = forecast2vec(cBase.forecast[i], h)
        cBaseDat[counter] = cBase.date[i]
        cBaseMod[counter] = cBase.model[i]
        cBaseLoc[counter] = cBase.location[i]
        cBaseHor[counter] = h
        counter += 1
    end
end

cBaseDF = [DataFrame(date = cBaseDat, model = cBaseMod, location = cBaseLoc, horizon = cBaseHor) DataFrame(cBaseMat, :auto)]
# CSV.write("cBase.csv", cBaseDF)