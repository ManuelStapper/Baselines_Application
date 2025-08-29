#############################################
### Creating forecasts of baseline models ###
#############################################

# Structure:
# 1) Load forecast data
# 2) Load baseline forecast data
# 3) Include some additional functions
# 4) Assess calibration
# 4a) How calibrated are the baselines?
# 4b) How uniformly calibrated are baselines?
# 4c) Is both even possible?
# 4d) How much does re-calibration help?
# 4e) Would it help to select different models for different locations? (maybe even horizons?)

# 5) Assess forecasts
# 5a) How does baseline selection affect the overall performance assessment?
# 5b) How does baseline selection affect the ranking of models?
# 5c) Contrast the above by two different performance measures. Which determines the ranking/overall performance more?
# 5d) Does the baseline change where a model is said to perform well vs poorly?
#     Select one model and few baselines with differnt calibratedness
# 5e) ...

using CSV, DataFrames, Dates, Plots
using BaselineModels
using StatsBase, Random, Optim, Distributions


# Set path to data here
path = "/Users/lshms101/Desktop/Projects/Baselines/Data"

# Load and prepare baseline forecasts

iBaseRaw = CSV.read(path * "/iBase.csv", DataFrame)
dates = unique(iBaseRaw.date)
locations = unique(iBaseRaw.location)
models = unique(iBaseRaw.model)

N = nrow(unique(iBaseRaw[:, 1:3]))

iBase = DataFrame(date = fill(Date(2000, 1, 1), N),
                  model = fill(1, N),
                  location = fill(1, N),
                  forecast = Vector{forecast}(undef, N))
#

dates1 = repeat(dates, length(locations))
locations1 = repeat(locations, inner = length(dates))

counter = 1
for i = models
    iBase.date[counter:counter + 4835] = dates1
    iBase.location[counter:counter + 4835] = locations1
    iBase.model[counter:counter + 4835] = fill(i, 4836)
    counter = counter + 4836
end


quantiles = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99]
α = quantiles[1:11]*2

function getBase(df, l, m, d)
    ind = (df.location .== l) .& (df.model .== m) .& (df.date .== d)
    x = df[ind, :]
    hSeq = x.horizon
    fi = (h -> forecastInterval(α, Vector{Float64}(x[h, 7:17]), Vector{Float64}(x[h, 18:end]))).(1:nrow(x))
    out = forecast(hSeq, median = x.x1, truth = x.x2, interval = fi)

    out
end

for i = 1:nrow(iBase)
    iBase.forecast[i] = getBase(iBaseRaw, iBase.location[i], iBase.model[i], iBase.date[i])
end

cBaseRaw = CSV.read(path * "/cBase.csv", DataFrame)
N = nrow(unique(cBaseRaw[:, 1:3]))
cBase = DataFrame(date = fill(Date(2000, 1, 1), N),
                  model = fill(1, N),
                  location = fill(1, N),
                  forecast = Vector{forecast}(undef, N))
#

dates = sort(unique(cBaseRaw.date))
models = sort(unique(cBaseRaw.model))
locations = sort(unique(cBaseRaw.location))

dates1 = repeat(dates, length(locations))
locations1 = repeat(locations, inner = length(dates))

counter = 1
for i = models
    cBase.date[counter:counter + 7903] = dates1
    cBase.location[counter:counter + 7903] = locations1
    cBase.model[counter:counter + 7903] = fill(i, 7904)
    counter = counter + 7904
end

for i = 1:nrow(cBase)
    cBase.forecast[i] = getBase(cBaseRaw, cBase.location[i], cBase.model[i], cBase.date[i])
end

# Load and prepare model forecasts
# Set path to code here
path = "/Users/lshms101/Desktop/Projects/Baselines/Code"

include(path*"/InputPrep.jl")
cd(path)

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

locations = sort(unique(ifc.location))

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

l = 11
d1 = covidTruth.date[covidTruth.location .== locations[l]]
t1 = covidTruth.value[covidTruth.location .== locations[l]]

covidTruth.value[(covidTruth.date .>= d1[112]) .& (covidTruth.location .== locations[11])] .= floor.(Int64, (t1[111:end-1] .+ t1[112:end])./2)

# Negative values in COVID truth data
# After visual checks, plausible that it is just a sign error
covidTruth.value = abs.(covidTruth.value)

# Check for positivity

function makePos(x::Vector{T}) where {T <: Real}
    if minimum(x) < 0
        x[x .< 0] .= 0
    end
    x
end

for i = 1:nrow(ifc)
    ifc.forecast[i].median = makePos(ifc.forecast[i].median)
    ifc.forecast[i].truth = abs.(ifc.forecast[i].truth)
    for h = 1:length(ifc.forecast[i].horizon)
        ifc.forecast[i].interval[h].l = makePos(ifc.forecast[i].interval[h].l)
        ifc.forecast[i].interval[h].u = makePos(ifc.forecast[i].interval[h].u)
    end
end

for i = 1:nrow(cfc)
    cfc.forecast[i].median = makePos(cfc.forecast[i].median)
    cfc.forecast[i].truth = abs.(cfc.forecast[i].truth)
    for h = 1:length(cfc.forecast[i].horizon)
        cfc.forecast[i].interval[h].l = makePos(cfc.forecast[i].interval[h].l)
        cfc.forecast[i].interval[h].u = makePos(cfc.forecast[i].interval[h].u)
    end
end

for i = 1:nrow(iBase)
    iBase.forecast[i].median = makePos(iBase.forecast[i].median)
    iBase.forecast[i].truth = abs.(iBase.forecast[i].truth)
    for h = 1:length(iBase.forecast[i].horizon)
        iBase.forecast[i].interval[h].l = makePos(iBase.forecast[i].interval[h].l)
        iBase.forecast[i].interval[h].u = makePos(iBase.forecast[i].interval[h].u)
    end
end

for i = 1:nrow(cBase)
    cBase.forecast[i].median = makePos(cBase.forecast[i].median)
    cBase.forecast[i].truth = abs.(cBase.forecast[i].truth)
    for h = 1:length(cBase.forecast[i].horizon)
        cBase.forecast[i].interval[h].l = makePos(cBase.forecast[i].interval[h].l)
        cBase.forecast[i].interval[h].u = makePos(cBase.forecast[i].interval[h].u)
    end
end


### Assess calibration

# Set path to general folder here
path = "/Users/lshms101/Desktop/Projects/Baselines/"

models = sort(unique(iBase.model))
modelsC = sort(unique(cBase.model))

modelNames = ["Constant", "Marginal", "KDE", "LSD", "OLS", "IDS",
              "ARMA", "INARCH", "ETS", "STL", "Constant (log)", "Marginal (log)",
              "KDE (log)", "LSD (log)", "OLS (log)", "IDS (log)", "ARMA (log)",
              "ETS (log)", "STL (log)"]
#

iPITfuns = Function[]
iDivs = Float64[]
cPITfuns = Function[]
cDivs = Float64[]

for i = 1:length(models)
    temp = iBase.forecast[(iBase.model .== i)]
    push!(iDivs, CvMdivergence(temp, [1, 2, 3, 4]))
    push!(iPITfuns, PITfun(temp, [1, 2, 3, 4]))
end

for i = 1:length(modelsC)
    temp = cBase.forecast[(cBase.model .== i)]
    push!(cDivs, CvMdivergence(temp, [1, 2, 3, 4]))
    push!(cPITfuns, PITfun(temp, [1, 2, 3, 4]))
end

iBest = argmin(iDivs)
iBest2 = findfirst(iDivs .== sort(iDivs)[2])
iWorst = argmax(iDivs)

cBest = argmin(cDivs)
cBest2 = findfirst(cDivs .== sort(cDivs)[2])
cWorst = argmax(cDivs)

p1 = plot(x -> x, xlim = (0, 1), label = "Perfect calibration", xlab = "u",
     ylab = "F(u)", color = "black", lw = 2, linestyle = :dash, legend = false)
for i = 1:length(iPITfuns)
    col = "gray"
    a = 0.3
    lab = ""
    if i == iBest
        col = "green"
        a = 1
        lab = modelNames[iBest]
    end
    if i == iWorst
        col = "red"
        a = 1
        lab = modelNames[iWorst]
    end
    plot!(iPITfuns[i], color = col, lw = 2, alpha = a, label = lab)
end
title!("Influenza")
annotate!(0.1, 0.9, text(modelNames[iBest], :green, 12, :left))
annotate!(0.1, 0.7, text(modelNames[iWorst], :red, 12, :left))
txt = "CvM = " * string(Int64(round(iDivs[iBest], digits = 0)))
annotate!(0.1, 0.85, text(txt, :green, :left, 12))
txt = "CvM = " * string(Int64(round(iDivs[iWorst], digits = 0)))
annotate!(0.1, 0.65, text(txt, :red, :left, 12))


p2 = plot(x -> x, xlim = (0, 1), label = "Perfect calibration", xlab = "u",
     ylab = "F(u)", color = "black", lw = 2, linestyle = :dash, legend = false)
for i = 1:length(cPITfuns)
    col = "gray"
    a = 0.3
    lab = ""
    if i == cBest
        col = "green"
        a = 1
        lab = modelNames[cBest]
    end
    if i == cWorst
        col = "red"
        a = 1
        lab = modelNames[cWorst]
    end
    plot!(cPITfuns[i], color = col, lw = 2, alpha = a, label = lab)
end
title!("COVID-19")
annotate!(0.65, 0.4, text(modelNames[cBest], :green, 12, :left))
annotate!(0.65, 0.2, text(modelNames[cWorst], :red, 12, :left))
txt = "CvM = " * string(Int64(round(cDivs[cBest], digits = 0)))
annotate!(0.65, 0.35, text(txt, :green, :left, 12))
txt = "CvM = " * string(Int64(round(cDivs[cWorst], digits = 0)))
annotate!(0.65, 0.15, text(txt, :red, :left, 12))

plot(p1, p2, layout = (1, 2), size = (900, 500), dpi = 300, left_margin = 5Plots.mm,
     bottom_margin = 5Plots.mm)
savefig(path*"Plots/PIT.png")
savefig(path*"Plots/PIT.pdf")


# Complete plot

plot(xlims = (-0.05, 5.05), ylims = (-0.05, 4.05), size = (1000, 800), axis = false, grid = false, legend = false, dpi = 300)
for i = 0:4
    plot!([0, 5], [i, i], color = "black", lw = ifelse(i == 2, 3, 1))
end
for i = 0:5
    plot!([i, i], [0, 4], color = "black")
end
for i = 1:5
    for j = 1:4
        plot!([i-1, i], [j-1, j], color = "black", linestyle = :dash)
    end
end
Plots.annotate!(-0.1, 3, Plots.text("Influenza", rotation = 90, 12))
Plots.annotate!(-0.1, 1, Plots.text("COVID-19", rotation = 90, 12))
xseq = collect(0:0.01:1)
xOffs = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
yOffsI = [3, 3, 3, 3, 3, 2, 2, 2, 2, 2]
yOffsC = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
for i = 1:10
    plot!(xseq .+ xOffs[i], iPITfuns[i].(xseq) .+ yOffsI[i], color = 1, lw = 2)
    txt = string(round(Int64, iDivs[i]))
    if i != 8
        plot!(xseq .+ xOffs[i], iPITfuns[i + 10 - (i > 8)].(xseq) .+ yOffsI[i], color = 1, lw = 2, linestyle = :dash)
        txt = txt*" ("*string(round(Int64, iDivs[i + 10 - (i > 8)]))*")"
    end
    Plots.annotate!(0.03 + xOffs[i], yOffsI[i] + 0.9, Plots.text(modelNames[i], :left))
    Plots.annotate!(0.03 + xOffs[i], yOffsI[i] + 0.8, Plots.text(txt, :left, 12))
end
for i = 1:10
    plot!(xseq .+ xOffs[i], cPITfuns[i].(xseq) .+ yOffsC[i], color = 1, lw = 2)
    txt = string(round(Int64, cDivs[i]))
    if i != 8
        plot!(xseq .+ xOffs[i], cPITfuns[i + 10 - (i > 8)].(xseq) .+ yOffsC[i], color = 1, lw = 2, linestyle = :dash)
        txt = txt*" ("*string(round(Int64, cDivs[i + 10 - (i > 8)]))*")"
    end
    Plots.annotate!(0.03 + xOffs[i], yOffsC[i] + 0.9, Plots.text(modelNames[i], :left))
    Plots.annotate!(0.03 + xOffs[i], yOffsC[i] + 0.8, Plots.text(txt, :left, 12))
end
savefig(path*"/Plots/PITall.png")
savefig(path*"/Plots/PITall.pdf")


### More detailled look into the PIT functions for best baselines
# Combine the plots into one!

p1 = plot(x -> x, xlim = (0, 1), label = "", xlab = "u", ylab = "F(u)", color = "black", lw = 2,
          linestyle = :dash,foreground_color_legend = nothing,
          background_color_legend = :transparent, left_margin = 10Plots.mm)
for h = 1:4
    temp = iBase.forecast[(iBase.model .== iBest)]
    plot!(PITfun(temp, [h]), label = "h = "*string(h), lw = 1)
end
title!("By Horizon")
annotate!(-0.4, 0.5, text("Influenza", 12, rotation = 90))

p2 = plot(x -> x, xlim = (0, 1), label = "", xlab = "u", ylab = "F(u)", color = "black", lw = 2,
          linestyle = :dash,foreground_color_legend = nothing,
          background_color_legend = :transparent)
for l = locations
    temp = iBase.forecast[(iBase.model .== iBest) .& (iBase.location .== l)]
    plot!(PITfun(temp, [1, 2, 3, 4]), label = "", lw = 1, color = "gray", alpha = 0.2)
end
title!("By Location")

p3 = plot(x -> x, xlim = (0, 1), label = "", xlab = "u", ylab = "F(u)", color = "black", lw = 2,
          linestyle = :dash,foreground_color_legend = nothing,
          background_color_legend = :transparent, left_margin = 10Plots.mm)
for h = 1:4
    temp = cBase.forecast[(cBase.model .== cBest)]
    plot!(PITfun(temp, [h]), label = "h = "*string(h), lw = 1)
end
annotate!(-0.4, 0.5, text("COVID-19", 12, rotation = 90))

p4 = plot(x -> x, xlim = (0, 1), label = "", xlab = "u", ylab = "F(u)", color = "black", lw = 2,
          linestyle = :dash,foreground_color_legend = nothing,
          background_color_legend = :transparent)
for l = locations
    temp = cBase.forecast[(cBase.model .== cBest) .& (cBase.location .== l)]
    plot!(PITfun(temp, [1, 2, 3, 4]), label = "", lw = 1, color = "gray", alpha = 0.2)
end

plot(p1, p2, p3, p4, layout = (2, 2), dpi = 300)
savefig(path*"/Plots/PITlh.png")
savefig(path*"/Plots/PITlh.pdf")



# What if we selected different models for each location
iBestVec = zeros(Int64, length(locations))
for l = 1:length(locations)
    iBestVec[l] = argmin((m -> CvMdivergence(iBase.forecast[(iBase.model .== m) .& (iBase.location .== locations[l])])).(models))
end

temp = vcat((l -> iBase.forecast[(iBase.model .== iBestVec[l]) .& (iBase.location .== locations[l])]).(1:length(locations))...)
plot(x -> x, xlim = (0, 1), label = "Perfect calibration", xlab = "u",
     ylab = "F(u)", color = "black", lw = 2, linestyle = :dash, dpi = 300)
plot!(PITfun(iBase.forecast[(iBase.model .== iBest)], [1, 2, 3, 4]), lw = 2, label = "Global best")
plot!(PITfun(temp, [1, 2, 3, 4]), lw = 2, label = "Individual best")

cBestVec = zeros(Int64, length(locations))
for l = 1:length(locations)
    cBestVec[l] = argmin((m -> CvMdivergence(cBase.forecast[(cBase.model .== m) .& (cBase.location .== locations[l])])).(models))
end
temp = vcat((l -> cBase.forecast[(cBase.model .== cBestVec[l]) .& (cBase.location .== locations[l])]).(1:length(locations))...)
plot(x -> x, xlim = (0, 1), label = "Perfect calibration", xlab = "u",
     ylab = "F(u)", color = "black", lw = 2, linestyle = :dash, dpi = 300)
plot!(PITfun(cBase.forecast[(cBase.model .== cBest)], [1, 2, 3, 4]), lw = 2, label = "Global best")
plot!(PITfun(temp, [1, 2, 3, 4]), lw = 2, label = "Individual best")

# Adding the region optimised baselines to collection
# ID = 20
for l = 1:length(locations)
    ind = (iBase.location .== locations[l]) .& (iBase.model .== iBestVec[l])
    dfAdd = DataFrame(date = iBase.date[ind], model = fill(20, 93), location = fill(locations[l], 93), forecast = iBase.forecast[ind])
    iBase = [iBase; dfAdd]    
end

for l = 1:length(locations)
    ind = (cBase.location .== locations[l]) .& (cBase.model .== cBestVec[l])
    dfAdd = DataFrame(date = cBase.date[ind], model = fill(20, sum(ind)), location = fill(locations[l], sum(ind)), forecast = cBase.forecast[ind])
    cBase = [cBase; dfAdd]   
end


# cols = [cols[1:10, 1]; cols[1:7, 2]; cols[9:10, 2]]
bestVecUni = sort(unique([iBestVec; cBestVec]))

cols = fill(RGB(1, 1, 1), 19)
cols[bestVecUni] = [
    RGB(170/255, 110/255, 40/255),   # Brown
    RGB(100/255, 100/255, 100/255),   # Yellow
    RGB(70/255, 240/255, 240/255),   # Cyan
    RGB(245/255, 130/255, 48/255),   # Orange
    RGB(60/255, 180/255, 75/255),    # Green
    RGB(145/255, 30/255, 180/255),   # Purple
    RGB(230/255, 25/255, 75/255),    # Red
    RGB(0/255, 130/255, 200/255),    # Blue X 
    RGB(128/255, 128/255, 0/255)     # Olive 
]
#

# Reading in Hex Tiles
# using Plots, Shapefile, GeoDataFrames, CSV, DataFrames, LibGEOS, StatsBase, GeoInterface, ColorSchemes
# shpHex = Shapefile.Table("/Users/lshms101/Desktop/Projects/Baselines/Data/Shapefile/StatesHex.shp")


colsVec = cols[iBestVec]
p1 = plot(shpHex.geometry[1], color = colsVec[1], aspect_ratio = 1, grid = false,
          xaxis = false, yaxis = false)
for i = 2:52
    plot!(shpHex.geometry[i], color = colsVec[i])
end
annotate!(shpHex.xCen, shpHex.yCen, Plots.text.(shpHex.abbr, color = "white", 12))

for i = 1:9
    annotate!([2, 4, 8][ceil(Int64, i / 3)], [5.5, 6.5, 6][(i % 3) + 1],
              Plots.text(modelNames[bestVecUni[i]], color = cols[bestVecUni[i]], :left, 12))
end

annotate!(0, 3, text("Influenza", 12, rotation = 90))

colsVec = cols[cBestVec]
p2 = plot(shpHex.geometry[1], color = colsVec[1], aspect_ratio = 1, grid = false,
          xaxis = false, yaxis = false)
for i = 2:52
    plot!(shpHex.geometry[i], color = colsVec[i])
end
annotate!(shpHex.xCen, shpHex.yCen, Plots.text.(shpHex.abbr, color = "white", 12))
annotate!(0, 3, text("COVID-19", 12, rotation = 90))

plot(p1, p2, layout = (2, 1), size = 0.7 .* (800, 1000), dpi = 300)

savefig(path*"/Plots/BestBaseline.png")
savefig(path*"/Plots/BestBaseline.pdf")

# Next also add the marginal model that uses future data
# only from times where forecasts are available from hubs

function horizon(x::forecast)
    x.horizon[end]
end

datesFci = unique(ifc.date)
hsi = (d -> maximum(horizon.(ifc.forecast[ifc.date .== d]))).(datesFci)
datesTrainingi = sort(unique(vcat((i -> datesFci[i] + Day(7).*collect(1:hsi[i])).(1:length(hsi))...)))

datesFcc = unique(cfc.date)
hsc = (d -> maximum(horizon.(cfc.forecast[cfc.date .== d]))).(datesFcc)
datesTrainingc = sort(unique(vcat((i -> datesFcc[i] + Day(7).*collect(1:hsc[i])).(1:length(hsc))...)))

function makeForecast(truth::DataFrame,
                      date::Date,
                      location::Int64,
                      datesTraining::Vector{Date},
                      hMax::Int64)
    #
    indL = (truth.location .== location)
    indTraining = (d -> d in datesTraining).(truth.date)
    ind = indL .& indTraining
    temp = truth[ind, :]

    meanVec = zeros(hMax)
    medianVec = zeros(hMax)
    intervalVec = Vector{forecastInterval}(undef, hMax)
    truthVec = zeros(hMax)
    for hh = 1:hMax
        targetDate = date + Day(7)*hh
        xx = temp.value[temp.date .!= targetDate]
        meanVec[hh] = median(xx)
        medianVec[hh] = median(xx)
        qVec = quantile(xx, quantiles)
        intervalVec[hh] = forecastInterval(α, qVec[1:11], qVec[end:-1:end-10])
        if any(temp.date .== targetDate)
            truthVec[hh] = temp.value[temp.date .== targetDate][1]
        else
            return forecast(collect(1:hh-1), mean = meanVec[1:hh-1], median = medianVec[1:hh-1], interval = intervalVec[1:hh-1], truth = truthVec[1:hh-1])        
        end
    end
    return forecast(collect(1:hMax), mean = meanVec, median = medianVec, interval = intervalVec, truth = truthVec)
end

for i = 1:length(datesFci)
    for l = 1:length(locations)
        fc = makeForecast(influenzaTruth, datesFci[i], locations[l], datesTrainingi, hsi[i])
        dfAdd = DataFrame(date = datesFci[i], model = 21, location = locations[l], forecast = fc)
        iBase = [iBase; dfAdd]
    end
end

for i = 1:length(datesFcc)
    for l = 1:length(locations)
        fc = makeForecast(covidTruth, datesFcc[i], locations[l], datesTrainingc, hsc[i])
        dfAdd = DataFrame(date = datesFcc[i], model = 21, location = locations[l], forecast = fc)
        cBase = [cBase; dfAdd]
    end
end

using StatsPlots
calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 17) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p11 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "", ylab = "F(u)", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
              markersize = 3, markerstrokewidth = 0)
end
title!("Best Single Baseline")

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 20) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p12 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "", ylab = "", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
    markersize = 3, markerstrokewidth = 0)
end
title!("Individual Best Baseline")

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 21) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p13 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "", ylab = "", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
    markersize = 3, markerstrokewidth = 0)
end
title!("Marginal (all)")

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 11) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p21 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "u", ylab = "F(u)", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
    markersize = 3, markerstrokewidth = 0)
end

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 20) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p22 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "u", ylab = "", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
    markersize = 3, markerstrokewidth = 0)
end

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 21) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

p23 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash, xlab = "u", ylab = "", lw = 2)
for i = 2:length(quantiles)
    boxplot!([quantiles[i]], calMatL[:, i], color = 1, bar_width = 0.03, alpha = 0.5,
    markersize = 3, markerstrokewidth = 0)
end
plot(p11, p12, p13, p21, p22, p23, layout=(2,3), size = (1000, 600), dpi = 300)
# savefig(path*"/Plots/UniformCalibration.png")


# Alternative with shapes instead of boxplots

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 17) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end
l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p11 = plot(x -> x, xlim = (0, 1), label = "", color = "black", linestyle = :dash,
           xlab = "", ylab = "F(u)", lw = 2, left_margin = 10Plots.mm)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "Range")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([3, 4, 4, 3], [0, 0, 1, 1]), color = 1, alpha = 0.6, linewidth = 0,
      label = "80% Interval")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "Median")
title!("Best Single Baseline")
annotate!(-0.5, 0.5, text("Influenza", 12, rotation = 90))


calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 20) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p12 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash,
           xlab = "", ylab = "", lw = 2)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "")
title!("Best Individual Baseline")

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(iBase.forecast[(iBase.model .== 21) .& (iBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p13 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash,
           xlab = "", ylab = "", lw = 2)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "")
title!("Marginal (all)")



calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 11) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p21 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash,
           xlab = "u", ylab = "F(u)", lw = 2)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "")
annotate!(-0.5, 0.5, text("COVID-19", 12, rotation = 90))

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 20) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p22 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash,
           xlab = "u", ylab = "", lw = 2)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "")

calMatL = zeros(length(locations), length(quantiles))
for l = 1:length(locations)
    calMatL[l, :] = PITfun(cBase.forecast[(cBase.model .== 21) .& (cBase.location .== locations[l])], [1, 2, 3, 4]).(quantiles)
end

l = [0; (i -> quantile(calMatL[:, i], 0.1)).(1:23); 1]
u = [0; (i -> quantile(calMatL[:, i], 0.9)).(1:23); 1]

l2 = [0; minimum(calMatL, dims = 1)[1, :]; 1]
u2 = [0; maximum(calMatL, dims = 1)[1, :]; 1]

p23 = plot(x -> x, xlim = (0, 1), legend = false, color = "black", linestyle = :dash,
           xlab = "u", ylab = "", lw = 2)
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l2; reverse(u2)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
plot!(Shape([0; quantiles; 1; 1; reverse(quantiles); 0], [l; reverse(u)]), color = 1, alpha = 0.3,
      linewidth = 0, label = "")
xlims!((0, 1))
plot!(quantiles, median(calMatL, dims = 1)[1, :], lw = 2, color = 2, label = "")

plot(p11, p12, p13, p21, p22, p23, layout=(2,3), size = 0.9 .* (1000, 600), left_margin = 12Plots.mm,
     bottom_margin = 5Plots.mm, dpi = 300)
#

savefig("/Users/lshms101/Desktop/Projects/Baselines/Plots/UniformCalibration2.png")
savefig("/Users/lshms101/Desktop/Projects/Baselines/Plots/UniformCalibration2.pdf")



#####################
### Testing block ###
#####################

# Tests for:
# 1) Exponential growth phases / log-trafo vs not
# 2) Seasonality
# 3) Global calibration
# 4) Uniform calibration (w.r.t. partitioning)


# 1) Log-transformation
using BoxCox

ibcs = (l -> fit(BoxCoxTransformation, influenzaTruth.value[influenzaTruth.location .== l] .+ 1).λ).(unique(influenzaTruth.location))
p1 = histogram(ibcs, title = "Influenza", label = "", xlabel = "λ")

cbcs = (l -> fit(BoxCoxTransformation, covidTruth.value[covidTruth.location .== l] .+ 1).λ).(unique(covidTruth.location))
p2 = histogram(cbcs, title = "COVID-19", label = "", xlabel = "λ")

plot(p1, p2)

iBCci = (l -> confint(fit(BoxCoxTransformation, influenzaTruth.value[influenzaTruth.location .== l] .+ 1))).(locations)
cBCci = (l -> confint(fit(BoxCoxTransformation, covidTruth.value[covidTruth.location .== l] .+ 1))).(locations)

iBCci = reshape(vcat(iBCci...), 2, 52)
cBCci = reshape(vcat(cBCci...), 2, 52)

sum(iBCci[1, :] .<= 0 .<= iBCci[2, :]) # For 12/52 regions log-transformation is in CI
sum(cBCci[1, :] .<= 0 .<= cBCci[2, :]) # For 0/52 regions log-transformation is in CI

plot(-1:0.01:1, (l -> sum(iBCci[1, :] .<= l .<= iBCci[2, :])).(-1:0.01:1), lw = 2, label = "Influenza")
plot!(-1:0.01:1, (l -> sum(cBCci[1, :] .<= l .<= cBCci[2, :])).(-1:0.01:1), lw = 2, label  = "COVID-19")


# 2) Seasonality

using LinearAlgebra

iACF = zeros(52, 101)
cACF = zeros(52, 101)
for i = 1:length(locations)
    x = influenzaTruth.value[influenzaTruth.location .== locations[i]]
    y = log.(x .+ 1)
    iACF[i, :] = (h -> cor(y[1 + h:end], y[1:end - h])).(0:100)
    
    x = covidTruth.value[covidTruth.location .== locations[i]]
    y = log.(x .+ 1)
    cACF[i, :] = (h -> cor(y[1 + h:end], y[1:end - h])).(0:100)
end

p1 = plot(0:100, mean(iACF, dims = 1)[1, :], lw = 2, legend = false, xlabel = "Lag", ylabel = "ACF", title = "Influenza")
plot!(0:100, (c -> quantile(iACF[:, c], 0.05)).(1:101), color = 1, linestyle = :dash)
plot!(0:100, (c -> quantile(iACF[:, c], 0.95)).(1:101), color = 1, linestyle = :dash)

p2 = plot(0:100, mean(cACF, dims = 1)[1, :], lw = 2, legend = false, color = 2, xlabel = "Lag", ylabel = "ACF", title = "COVID-19")
plot!(0:100, (c -> quantile(cACF[:, c], 0.05)).(1:101), color = 2, linestyle = :dash)
plot!(0:100, (c -> quantile(cACF[:, c], 0.95)).(1:101), color = 2, linestyle = :dash)

plot(p1, p2)


function trendSeason(y, per)
    T = length(y)
    X = [ones(T-1) sin.((1:T-1) ./ per .* (2*π)) cos.((1:T-1) ./ per .* (2*π)) collect(1:T-1) y[1:end-1]]
    beta = inv(X'X)*X'y[2:end]
    SSR = sum((y[2:end] .- X*beta).^2)
    Sigma = var(y[2:end] .- X*beta)*inv(X'X)

    # Test (positive) Autocorrelation
    p1 = 2 .* ccdf(TDist(T - 4 - 1), abs(beta[5] ./ sqrt.(Sigma[5, 5])))

    # Testing trend
    p2 = 2 .* ccdf(TDist(T - 4 - 1), abs(beta[4] ./ sqrt.(Sigma[4, 4])))

    # Testing seasonality jointly
    X2 = [ones(T-1) collect(1:T-1)]
    beta2 = inv(X2'X2)*X2'y[2:end]
    SSRr = sum((y[2:end] .- X2*beta2).^2)

    Fstat = ((SSRr - SSR)/2)/ (SSR/(T - 4 - 1))
    # p-Value of F test for seasonality
    p3 = ccdf(FDist(2, T - 4 - 1), Fstat)
    return [p1, p2, p3]
end

ips = (l -> trendSeason(log.(influenzaTruth.value[influenzaTruth.location .== l][45:end] .+ 1), 52)).(unique(influenzaTruth.location))
sum(first.(ips) .< 0.05)
sum((xx -> xx[2]).(ips) .< 0.05)
sum(last.(ips) .< 0.05)

cps = (l -> trendSeason(log.(covidTruth.value[covidTruth.location .== l] .+ 1), 52)).(unique(covidTruth.location))
sum(first.(cps) .< 0.05)
sum((xx -> xx[2]).(cps) .< 0.05)
sum(last.(cps) .< 0.05)


# 3) Global Calibration

# What is the distribution of CvM under H₀: Perfect calibration?
# Example: One baseline
b = 1
locVec = iBase.location[iBase.model .== b]
fcs = iBase.forecast[iBase.model .== b]
hVec = (x -> x.horizon).(fcs)
hs = vcat(hVec...)
locs = vcat(fill.(locVec, length.(hVec))...)
steps = vcat((i -> (hh -> makeStep(fcs[i], hh)).(hVec[i])).(1:length(fcs))...)

# Depends on number of forecasts (maybe there is a limiting distribution?)

p = round.(diff([0.0; quantiles; 1.0]), digits = 4)
0.416 - sum(p.^3)/12 # 95\% quantile

tStat = zeros(10000)

sam = rand(Multinomial(length(steps), p))

function CvMdivergence2(sam::Vector{Int64}, quantiles::Vector{Float64})
    αs = [0.0; quantiles; 1.0]
    τs = [0.0; cumsum(sam) ./ sum(sam)]

    out = 0.0
    for i in 1:length(αs)-1
        α₁, α₂ = αs[i], αs[i+1]
        τ₁, τ₂ = τs[i], τs[i+1]

        h = α₂ - α₁
        m = (τ₂ - τ₁) / h
        b = τ₁ - m * α₁
        a = m - 1
        I = (a^2 / 3) * (α₂^3 - α₁^3) + (a * b) * (α₂^2 - α₁^2) + (b^2) * (α₂ - α₁)

        out += I
    end

    return sum(sam)*out
end

tStat = zeros(1000000)

for i = 1:1000000
    Random.seed!(i)
    sam = rand(Multinomial(length(steps), p))
    tStat[i] = CvMdivergence2(sam, quantiles)
end

histogram(tStat)
vline!([0.416 - sum(p.^3)/12])

quantile(tStat, 0.95)
quantile(tStat, 0.99)
quantile(tStat, 0.999)
quantile(tStat, 0.9999)

0.416 - sum(p.^3)/12 # Confirms integrated squared brownian bridge distribution

sort(iDivs)
sort(cDivs)

# 4) Uniform convergence
# With respect to a partitioning (example here: location)
# Two different test procedures: H₀: Uniform calibration or H₀: No difference between partitions and full sample

# Version 1:
# 52 locations of 366 forecasts each (influenza)
# Important for limiting distribution only the number of locations

iDivsLoc = zeros(21)
cDivsLoc = zeros(21)

for b = 1:21
    iDivsLoc[b] = sum((l -> CvMdivergence(iBase.forecast[(iBase.model .== b) .& (iBase.location .== l)], [1, 2, 3, 4])).(locations))
    cDivsLoc[b] = sum((l -> CvMdivergence(cBase.forecast[(cBase.model .== b) .& (cBase.location .== l)], [1, 2, 3, 4])).(locations))
end

sam2 = zeros(100000)
d0 = Multinomial(366, p)
for i = 1:100000
    Random.seed!(i)
    for l = 1:52
        sam2[i] += CvMdivergence2(rand(d0), quantiles)
    end
end

histogram(sam2)
mean(sam2)
quantile(sam2, 0.95)

sort(iDivsLoc)
sort(cDivsLoc)

quantile(sam2, 0.95)
sort(iDivsLoc)
sort(cDivsLoc)

sam3 = zeros(10000, 21, 2)

for b = 1:21
    Random.seed!(b)
    locVec = iBase.location[iBase.model .== b]
    fcs = iBase.forecast[iBase.model .== b]
    hVec = (x -> x.horizon).(fcs)
    hs = vcat(hVec...)
    locs = vcat(fill.(locVec, length.(hVec))...)
    steps = vcat((i -> (hh -> makeStep(fcs[i], hh)).(hVec[i])).(1:length(fcs))...)
    lowers = (x -> x.l).(steps)
    p = (q -> mean(lowers .== q)).([0.0; quantiles])
    d0 = Multinomial(366, p)

    for i = 1:10000
        for j = 1:52
            sam3[i, b, 1] += CvMdivergence2(rand(d0), quantiles)
        end
    end

    locVec = cBase.location[cBase.model .== b]
    fcs = cBase.forecast[cBase.model .== b]
    hVec = (x -> x.horizon).(fcs)
    hs = vcat(hVec...)
    locs = vcat(fill.(locVec, length.(hVec))...)
    steps = vcat((i -> (hh -> makeStep(fcs[i], hh)).(hVec[i])).(1:length(fcs))...)
    lowers = (x -> x.l).(steps)
    p = (q -> mean(lowers .== q)).([0.0; quantiles])
    d0 = Multinomial(602, p)

    for i = 1:10000
        for j = 1:52
            sam3[i, b, 2] += CvMdivergence2(rand(d0), quantiles)
        end
    end
end

scatter((b -> mean(sam3[:, b, 1] .> iDivsLoc[b])).(1:21))
hline!([0.05])
scatter((b -> mean(sam3[:, b, 2] .> cDivsLoc[b])).(1:21))
hline!([0.05])

(b -> quantile(sam3[:, b, 1], 0.025) .<= iDivsLoc[b] .<= quantile(sam3[:, b, 1], 0.975)).(1:21) |> findall

(b -> quantile(sam3[:, b, 2], 0.025) .<= cDivsLoc[b] .<= quantile(sam3[:, b, 2], 0.975)).(1:21) |> findall
[modelNames; "Indiv"; "Marginal (all)"][[6, 10, 16, 17, 19]]


b = 1
histogram(sam3[:, b, 1])
vline!([iDivsLoc[b]], lw = 2)

histogram(sam3[:, b, 2])
vline!([cDivsLoc[b]], lw = 2)
cDivs[b]

#########################
### Testing block end ###
#########################


### Evaluating forecast models
# Comparison of models with best vs worst baseline

function findA(x::DataFrame, y::DataFrame)
    nx = nrow(x)
    ny = nrow(y)
    
    keepx = fill(false, nrow(x))
    keepy = fill(false, nrow(y))

    for i = 1:nx
        ind = (y.date .== x.date[i]) .& (y.location .== x.location[i])
        if sum(ind) == 1
            keepx[i] = true
            keepy[findfirst(ind)] = true
        end
    end
    return x[keepx, :], y[keepy, :]
end

function relativeSkill(dati::DataFrame, datj::DataFrame, logTrafo::Bool = false)
    dati2 = deepcopy(dati)
    datj2 = deepcopy(datj)
    dati2 = sort(dati2, [1, 3])
    datj2 = sort(datj2, [1, 3])
    hi = getproperty.(dati2.forecast, :horizon)
    hj = getproperty.(datj2.forecast, :horizon)
    hs = intersect.(hi, hj)

    scorei = 0.0
    scorej = 0.0
    for k = 1:nrow(dati)
        scorei += sum((h -> WIS(dati2.forecast[k], h, logTrafo = logTrafo)).(hs[k]))
        scorej += sum((h -> WIS(datj2.forecast[k], h, logTrafo = logTrafo)).(hs[k]))
    end
    scorei/scorej
end

m = 1 # selecting the baseline
ifcRel = copy(ifc)
tempAdd = iBase[(iBase.model .== m), :]
tempAdd.model .= 0
ifcRel = [ifcRel; tempAdd]

M = length(unique(ifcRel.model))
relWIS = zeros(M, length(models), 2)

θ = ones(M, M, 2)
modelsRel = sort(unique(ifcRel.model))
for i = 1:M
    for j = i:M
        dati = ifcRel[ifcRel.model .== modelsRel[i], :]
        datj = ifcRel[ifcRel.model .== modelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[i, j, 1] = relativeSkill(dati, datj, false)
            θ[i, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, i, 1] = 1/θ[i, j, 1]
        θ[j, i, 2] = 1/θ[i, j, 2]
    end
end

θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
θ2 = θ2 ./ θ2[1]
relWIS[:, m, 1] = θ2

θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
θ2 = θ2 ./ θ2[1]
relWIS[:, m, 2] = θ2

models = collect(1:19)

for m = 2:length(models)
    ifcRel = copy(ifc)
    tempAdd = iBase[(iBase.model .== m), :]
    tempAdd.model .= 0
    ifcRel = [ifcRel; tempAdd]

    M = length(unique(ifcRel.model))
    modelsRel = sort(unique(ifcRel.model))
    for j = 1:M
        dati = ifcRel[ifcRel.model .== modelsRel[1], :]
        datj = ifcRel[ifcRel.model .== modelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[1, j, 1] = relativeSkill(dati, datj, false)
            θ[1, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, 1, 1] = 1/θ[1, j, 1]
        θ[j, 1, 2] = 1/θ[1, j, 2]
    end
    θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    relWIS[:, m, 1] = θ2

    θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    relWIS[:, m, 2] = θ2

    println(now())
    println(m)
end

relWISadd = zeros(82, 2, 2)

for m = 20:21
    ifcRel = copy(ifc)
    tempAdd = iBase[(iBase.model .== m), :]
    tempAdd.model .= 0
    ifcRel = [ifcRel; tempAdd]

    M = length(unique(ifcRel.model))
    modelsRel = sort(unique(ifcRel.model))
    for j = 1:M
        dati = ifcRel[ifcRel.model .== modelsRel[1], :]
        datj = ifcRel[ifcRel.model .== modelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[1, j, 1] = relativeSkill(dati, datj, false)
            θ[1, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, 1, 1] = 1/θ[1, j, 1]
        θ[j, 1, 2] = 1/θ[1, j, 2]
    end
    θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    relWISadd[:, m-19, 1] = θ2

    θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    relWISadd[:, m-19, 2] = θ2

    println(now())
    println(m)
end




m = 1 # selecting the baseline
cfcRel = copy(cfc)
tempAdd = cBase[(cBase.model .== m), :]
tempAdd.model .= 0
cfcRel = [cfcRel; tempAdd]

M = length(unique(cfcRel.model))
crelWIS = zeros(M, length(models), 2)

θ = ones(M, M, 2)
cmodelsRel = sort(unique(cfcRel.model))
for i = 1:M
    for j = i:M
        dati = cfcRel[cfcRel.model .== cmodelsRel[i], :]
        datj = cfcRel[cfcRel.model .== cmodelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[i, j, 1] = relativeSkill(dati, datj, false)
            θ[i, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, i, 1] = 1/θ[i, j, 1]
        θ[j, i, 2] = 1/θ[i, j, 2]
    end
end

θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
θ2 = θ2 ./ θ2[1]
crelWIS[:, m, 1] = θ2

θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
θ2 = θ2 ./ θ2[1]
crelWIS[:, m, 2] = θ2

for m = 2:length(models)
    cfcRel = copy(cfc)
    tempAdd = cBase[(cBase.model .== m), :]
    tempAdd.model .= 0
    cfcRel = [cfcRel; tempAdd]

    M = length(unique(cfcRel.model))
    cmodelsRel = sort(unique(cfcRel.model))
    for j = 1:M
        dati = cfcRel[cfcRel.model .== cmodelsRel[1], :]
        datj = cfcRel[cfcRel.model .== cmodelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[1, j, 1] = relativeSkill(dati, datj, false)
            θ[1, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, 1, 1] = 1/θ[1, j, 1]
        θ[j, 1, 2] = 1/θ[1, j, 2]
    end
    θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    crelWIS[:, m, 1] = θ2

    θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    crelWIS[:, m, 2] = θ2

    println(now())
    println(m)
end

crelWISadd = zeros(58, 2, 2)

for m = 20:21
    cfcRel = copy(cfc)
    tempAdd = cBase[(cBase.model .== m), :]
    tempAdd.model .= 0
    cfcRel = [cfcRel; tempAdd]

    M = length(unique(cfcRel.model))
    cmodelsRel = sort(unique(cfcRel.model))
    for j = 1:M
        dati = cfcRel[cfcRel.model .== cmodelsRel[1], :]
        datj = cfcRel[cfcRel.model .== cmodelsRel[j], :]
        dati, datj = findA(dati, datj)
        if nrow(dati) > 0
            θ[1, j, 1] = relativeSkill(dati, datj, false)
            θ[1, j, 2] = relativeSkill(dati, datj, true)
        end
        θ[j, 1, 1] = 1/θ[1, j, 1]
        θ[j, 1, 2] = 1/θ[1, j, 2]
    end
    θ2 = (i -> prod(θ[i, :, 1])^(1/(sum(θ[i, :, 1] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    crelWISadd[:, m-19, 1] = θ2

    θ2 = (i -> prod(θ[i, :, 2])^(1/(sum(θ[i, :, 2] .!= 1) + 1))).(1:size(θ)[1])
    θ2 = θ2 ./ θ2[1]
    crelWISadd[:, m-19, 2] = θ2

    println(now())
    println(m)
end

relWIS = [relWIS;; relWISadd]
crelWIS = [crelWIS;; crelWISadd]

# Save relWIS and crelWIS
# CSV.write("relWIS.csv", DataFrame([relWIS[:, :, 1]; relWIS[:, :, 2]], :auto))
# CSV.write("crelWIS.csv", DataFrame([crelWIS[:, :, 1]; crelWIS[:, :, 2]], :auto))

temp = CSV.read("relWIS.csv", DataFrame)
relWIS = zeros(82, 21, 2)
relWIS[:, :, 1] = Matrix(temp[1:82, :])
relWIS[:, :, 2] = Matrix(temp[83:end, :])

temp = CSV.read("crelWIS.csv", DataFrame)
crelWIS = zeros(58, 21, 2)
crelWIS[:, :, 1] = Matrix(temp[1:58, :])
crelWIS[:, :, 2] = Matrix(temp[59:end, :])

# Plots:
# Show absolute differences in relative skill due to baseline selection
# Show that ranking does not really change with baseline selection
# Show that ranking DOES change with evaluation metric
# With one metric: Depends on where the overlap is:
## If a baseline is great at "constant" times and one model mainly forecasts there, it's viewed less good.

# Influenza - original WIS
nBins = 20

temp = vec(relWIS[2:end, [argmin(iDivs), argmax(iDivs)], 1])
m = minimum(temp)
M = maximum(temp)
binSeq = collect(m:(M - m)/nBins:M)
p1 = histogram(relWIS[2:end, argmin(iDivs), 1], alpha = 0.5, bins = binSeq,
               label = modelNames[argmin(iDivs)], color = "green", legend = false)
histogram!(relWIS[2:end, argmax(iDivs), 1], alpha = 0.5, bins = binSeq,
           label = modelNames[argmax(iDivs)], color = "red")
vline!([1], lw = 2, color = "black", linestyle = :dash, label = "")
title!("WIS")
ylabel!("Influenza")
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.6,
          text(modelNames[argmin(iDivs)], :green, 8, :left))
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.4,
          text(modelNames[argmax(iDivs)], :red, 8, :left))

# Influenza WIS on logs
temp = vec(relWIS[2:end, [argmin(iDivs), argmax(iDivs)], 2])
m = minimum(temp)
M = maximum(temp)
binSeq = collect(m:(M - m)/nBins:M)
p2 = histogram(relWIS[2:end, argmin(iDivs), 2], alpha = 0.5, bins = binSeq,
               label = modelNames[argmin(iDivs)], color = "green", legend = false)
histogram!(relWIS[2:end, argmax(iDivs), 2], alpha = 0.5,  bins = binSeq,
           label = modelNames[argmax(iDivs)], color = "red")
vline!([1], lw = 2, color = "black", linestyle = :dash, label = "")
title!("lWIS")
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.6,
          text(modelNames[argmin(iDivs)], :green, 8, :left))
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.4,
          text(modelNames[argmax(iDivs)], :red, 8, :left))


# COVID - original WIS
temp = vec(crelWIS[2:end, [argmin(cDivs), argmax(cDivs)], 1])
m = minimum(temp)
M = maximum(temp)
binSeq = collect(m:(M - m)/nBins:M)
p3 = histogram(crelWIS[2:end, argmin(cDivs), 1], alpha = 0.5, bins = binSeq,
               label = modelNames[argmin(cDivs)], color = "green", legend = false)
histogram!(crelWIS[2:end, argmax(cDivs), 1], alpha = 0.5, bins = binSeq,
           label = modelNames[argmax(cDivs)], color = "red")
vline!([1], lw = 2, color = "black", linestyle = :dash, label = "")
title!("")
xlabel!("Scaled Relative Skill Score")
ylabel!("COVID-19")
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.6,
          text(modelNames[argmin(cDivs)], :green, 8, :left))
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.4,
          text(modelNames[argmax(cDivs)], :red, 8, :left))

# COVID - WIS on logs
temp = vec(crelWIS[2:end, [argmin(cDivs), argmax(cDivs)], 2])
m = minimum(temp)
M = maximum(temp)
binSeq = collect(m:(M - m)/nBins:M)
p4 = histogram(crelWIS[2:end, argmin(cDivs), 2], alpha = 0.5, bins = binSeq,
               label = modelNames[argmin(cDivs)], color = "green", legend = false)
histogram!(crelWIS[2:end, argmax(cDivs), 2], alpha = 0.5, bins = binSeq,
           label = modelNames[argmax(cDivs)], color = "red")
vline!([1], lw = 2, color = "black", linestyle = :dash, label = "")
title!("")
xlabel!("Scaled Relative Skill Score")
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.6,
          text(modelNames[argmin(cDivs)], :green, 8, :left))
annotate!(m + 0.62*(M - m), ylims(current())[2] * 0.4,
          text(modelNames[argmax(cDivs)], :red, 8, :left))

plot(p1, p2, p3, p4, layout=(2, 2), dpi = 300)

savefig(path*"/Plots/RelativeSkill1.png")
savefig(path*"/Plots/RelativeSkill1.pdf")


sum(relWIS[2:end, argmin(iDivs), 1] .< 1)
sum(relWIS[2:end, argmax(iDivs), 1] .< 1)

sum(relWIS[2:end, argmin(iDivs), 2] .< 1)
sum(relWIS[2:end, argmax(iDivs), 2] .< 1)

sum(crelWIS[2:end, argmin(cDivs), 1] .< 1)
sum(crelWIS[2:end, argmax(cDivs), 1] .< 1)

sum(crelWIS[2:end, argmin(cDivs), 2] .< 1)
sum(crelWIS[2:end, argmax(cDivs), 2] .< 1)


rMat = zeros(Int64, 81, 21, 2)
for m = 1:21
    temp = relWIS[2:end, m, 1]
    rMat[:, m, 1] = (i -> findfirst(temp[i] .== sort(temp))).(1:81)
    temp = relWIS[2:end, m, 2]
    rMat[:, m, 2] = (i -> findfirst(temp[i] .== sort(temp))).(1:81)
end

relWIS
crMat = zeros(Int64, 57, 21, 2)
for m = 1:21
    temp = crelWIS[2:end, m, 1]
    crMat[:, m, 1] = (i -> findfirst(temp[i] .== sort(temp))).(1:57)
    temp = crelWIS[2:end, m, 2]
    crMat[:, m, 2] = (i -> findfirst(temp[i] .== sort(temp))).(1:57)
end

modelNamesShort = copy(modelNames)
modelNamesShort[1] = "Const"
modelNamesShort[11] = "Const-l"
modelNamesShort[2] = "Mgn"
modelNamesShort[12] = "Mgn-l"

modelNamesShort[13] = "KDE-l"
modelNamesShort[14] = "LSD-l"
modelNamesShort[15] = "OLS-l"
modelNamesShort[16] = "IDS-l"
modelNamesShort[17] = "ARMA-l"
modelNamesShort[18] = "ETS-l"
modelNamesShort[19] = "STL-l"

cols = palette(:viridis, 81)

o = (i -> findfirst(iDivs .== sort(iDivs)[i])).(1:19)

p1 = plot(83 .- rMat[1, o, 1], legend = false, lw = 2, color = cols[Int64.(rMat[1, 17, 1])],
          left_margin = 15Plots.mm)
for i = 2:81
    plot!(83 .- rMat[i, o, 1], color = cols[Int64.(rMat[i, 17, 1])], lw = 2)
end
xlabel!("")
ylabel!("Rank of Forecast Model")
yticks!(83 .- collect(10:10:80), string.(10:10:80))
xticks!(1:19, modelNamesShort[(i -> findfirst(iDivs .== sort(iDivs)[i])).(1:19)], rotation = 45)
title!("WIS")
annotate!(-4.5, 40, text("Influenza", 12, rotation = 90))

p2 = plot(83 .- rMat[1, o, 2], legend = false, lw = 2, color = cols[Int64.(rMat[1, 17, 2])])
for i = 2:81
    plot!(83 .- rMat[i, o, 2], color = cols[Int64.(rMat[i, 17, 2])], lw = 2)
end
xlabel!("")
title!("lWIS")
yticks!(83 .- collect(10:10:80), string.(10:10:80))
xticks!(1:19, modelNamesShort[(i -> findfirst(iDivs .== sort(iDivs)[i])).(1:19)], rotation = 45)

cols = palette(:viridis, 57)
o = (i -> findfirst(cDivs .== sort(cDivs)[i])).(1:19)

p3 = plot(59 .- crMat[1, o, 1], legend = false, lw = 2, color = cols[Int64.(crMat[1, 11, 1])],
          left_margin = 15Plots.mm)
for i = 2:57
    plot!(59 .- crMat[i, o, 1], color = cols[Int64.(crMat[i, 11, 1])], lw = 2)
end
ylabel!("")
ylabel!("Rank of Forecast Model")
yticks!(59 .- collect(10:10:50), string.(10:10:50))
xticks!(1:19, modelNamesShort[(i -> findfirst(cDivs .== sort(cDivs)[i])).(1:19)], rotation = 45)
annotate!(-4.5, 27, text("COVID-19", 12, rotation = 90))

p4 = plot(59 .- crMat[1, o, 2], legend = false, lw = 2, color = cols[Int64.(crMat[1, 11, 2])])
for i = 2:57
    plot!(59 .- crMat[i, o, 2], color = cols[Int64.(crMat[i, 11, 2])], lw = 2)
end
xlabel!("")
ylabel!("")
yticks!(59 .- collect(10:10:50), string.(10:10:50))
xticks!(1:19, modelNamesShort[(i -> findfirst(cDivs .== sort(cDivs)[i])).(1:19)], rotation = 45)

plot(p1, p2, p3, p4, layout = (2, 2), size = 0.9 .* (800, 650), dpi = 300)
savefig(path*"/Plots/RankPlot.png")
savefig(path*"/Plots/RankPlot.pdf")


# Also create a heatmap of Kendall's tau for pairwise comparisons
CI = zeros(19, 19)
CIl = zeros(19, 19)
CC = zeros(19, 19)
CCl = zeros(19, 19)

iDivs2 = [iDivs; CvMdivergence(iBase.forecast[iBase.model .== 20]); CvMdivergence(iBase.forecast[iBase.model .== 21])]
cDivs2 = [cDivs; CvMdivergence(cBase.forecast[cBase.model .== 20]); CvMdivergence(cBase.forecast[cBase.model .== 21])]

for i = 1:19
    for j = 1:19
        bi = findfirst(iDivs .== sort(iDivs)[i])
        bj = findfirst(iDivs .== sort(iDivs)[j])
        CI[i, j] = corkendall(rMat[:, bi, 1], rMat[:, bj, 1])
        CIl[i, j] = corkendall(rMat[:, bi, 2], rMat[:, bj, 2])

        bi = findfirst(cDivs .== sort(cDivs)[i])
        bj = findfirst(cDivs .== sort(cDivs)[j])
        CC[i, j] = corkendall(crMat[:, bi, 1], crMat[:, bj, 1])
        CCl[i, j] = corkendall(crMat[:, bi, 2], crMat[:, bj, 2])
    end
end

cLower = 0.9

p1 = heatmap(CI[end:-1:1, :], grid = false, xaxis = false, yaxis = false, clims = (cLower, 1),
             legend = false, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
for i = 1:19
    temp = modelNamesShort[findfirst(iDivs .== sort(iDivs)[i])]
    annotate!(0, 20 - i, text(temp, 8, :right))
    annotate!(i, 0, text(temp, 8, :right, rotation = 60))
end
title!("WIS")
plot!([0.5, 19.5], [-3, -3], label = "", arrow = :closed, color = "black")
plot!([-3, -3], [19.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(10, -4, text("Calibration", 8))
annotate!(2, -4, text("Good", 8, :left))
annotate!(18, -4, text("Poor", 8, :right))

annotate!(-6, 10, text("Influenza", 12, rotation = 90))
annotate!(-4.5, 10, text("Calibration", 8, rotation = 90))
annotate!(-4.5, 18, text("Good", 8, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 8, :left, rotation = 90))

p2 = heatmap(CIl[end:-1:1, :], grid = false, xaxis = false, yaxis = false, clims = (cLower, 1),
             legend = false, bottom_margin = 5Plots.mm)
for i = 1:19
    temp = modelNamesShort[findfirst(iDivs .== sort(iDivs)[i])]
    annotate!(0, 20 - i, text(temp, 8, :right))
    annotate!(i, 0, text(temp, 8, :right, rotation = 60))
end
title!("lWIS")
plot!([0.5, 19.5], [-3, -3], label = "", arrow = :closed, color = "black")
plot!([-3, -3], [19.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(10, -4, text("Calibration", 8))
annotate!(2, -4, text("Good", 8, :left))
annotate!(18, -4, text("Poor", 8, :right))

annotate!(-4.5, 10, text("Calibration", 8, rotation = 90))
annotate!(-4.5, 18, text("Good", 8, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 8, :left, rotation = 90))

p3 = heatmap(CC[end:-1:1, :], grid = false, xaxis = false, yaxis = false, clims = (cLower, 1),
             legend = false, left_margin = 10Plots.mm, bottom_margin = 5Plots.mm)
for i = 1:19
    temp = modelNamesShort[findfirst(cDivs .== sort(cDivs)[i])]
    annotate!(0, 20 - i, text(temp, 8, :right))
    annotate!(i, 0, text(temp, 8, :right, rotation = 60))
end
plot!([0.5, 19.5], [-3, -3], label = "", arrow = :closed, color = "black")
plot!([-3, -3], [19.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(10, -4, text("Calibration", 8))
annotate!(2, -4, text("Good", 8, :left))
annotate!(18, -4, text("Poor", 8, :right))

annotate!(-6, 10, text("COVID-19", 12, rotation = 90))
annotate!(-4.5, 10, text("Calibration", 8, rotation = 90))
annotate!(-4.5, 18, text("Good", 8, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 8, :left, rotation = 90))

p4 = heatmap(CCl[end:-1:1, :], grid = false, xaxis = false, yaxis = false, clims = (cLower, 1),
             legend = false, bottom_margin = 5Plots.mm)
for i = 1:19
    temp = modelNamesShort[findfirst(cDivs .== sort(cDivs)[i])]
    annotate!(0, 20 - i, text(temp, 8, :right))
    annotate!(i, 0, text(temp, 8, :right, rotation = 60))
end
plot!([0.5, 19.5], [-3, -3], label = "", arrow = :closed, color = "black")
plot!([-3, -3], [19.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(10, -4, text("Calibration", 8))
annotate!(2, -4, text("Good", 8, :left))
annotate!(18, -4, text("Poor", 8, :right))

annotate!(-4.5, 10, text("Calibration", 8, rotation = 90))
annotate!(-4.5, 18, text("Good", 8, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 8, :left, rotation = 90))

dummy_data = fill(missing, 2, 2)
p_colorbar = heatmap(dummy_data, clims = (cLower, 1.0),
 colorbar = :right,
 showaxis = false,
 grid = false,
 framestyle = :none,
 background_color = :transparent,
 foreground_color = :transparent,
 ytickfontsize = 10)
#
pEmpty = plot([1, 2], legend = false, axis = false, grid = false, color = :transparent)
l1 = @layout [a{0.1h}; b{0.8h}; c{0.1h}]
pBar = plot(pEmpty, p_colorbar, pEmpty, layout = l1, size = (100, 900))
# Use a custom layout
l = @layout [[a b; c d] e{0.05w}]
plot(p1, p2, p3, p4, pBar, layout = l, size = 0.9 .* (1000, 900), dpi = 300)
savefig(path*"/Plots/RankHeatmap.png")
savefig(path*"/Plots/RankHeatmap.pdf")


# Look at extreme correlation:
# COVID-19 best two baselines are different in their rankings

scatter(crMat[:, 17, 1], crMat[:, 15, 1])
scatter(crMat[:, 17, 1], crMat[:, 17, 1] .- crMat[:, 15, 1])
argmax((crMat[:, 17, 1] .- crMat[:, 15, 1]))
argmin((crMat[:, 17, 1] .- crMat[:, 15, 1]))

# Models 40 and 46 have the largest differences in rank between the two baselines
# One extreme in one direction and one in the other
# Using ARMA-l: Model 40 ranked 45th, using OLS-l 31st
# Model 40 only submitted forecasts during strong increase at the end of 2021/beginning of 2022
# OLS may overshoot (especially since fitted in logs), making model look better
# ARMA not as affected by sudden increase, hence more realistic comparison (?)
# Using ARMA-l: Model 46 ranked 30th, using OLS-l 35st
# Model 46 only submitted forecasts from the decline at the beginning of 2022
# ARMA-l maybe still affected by large numbers before (slowly returning to normal levels)
# OLS reacts more quickly?

# MUNI-VAR: Largest difference in ranks
crMat[40, 17, 1]
crMat[40, 15, 1]

# IHME-CurveFit: Second largest difference (same direction)
crMat[21, 17, 1]
crMat[21, 15, 1]

# SDSC_ISG-TrendModel (largest difference in the other direction)
crMat[46, 17, 1]
crMat[46, 15, 1]

h = 1

datesTemp = sort(unique(cfc.date[cfc.model .== 46]))
indTemp = findall((covidTruth.location .== 0) .& ((d -> d in datesTemp).(covidTruth.date)))
plot(covidTruth.date[covidTruth.location .== 0], covidTruth.value[covidTruth.location .== 0],
     label = "Truth", xlabel = "Date", ylabel = "Cases", lw = 3)
cfcTemp = cfc[(cfc.model .== 46) .& (cfc.location .== 0), :]
baseTemp1 = cBase[(cBase.model .== 17) .& (cBase.location .== 0), :]
baseTemp2 = cBase[(cBase.model .== 15) .& (cBase.location .== 0), :]
scatter!(cfcTemp.date .+ Day(7*h), (x -> x.median[h]).(cfcTemp.forecast), lw = 2,
      label = "Forecast model")
plot!(baseTemp1.date[1:end-h] .+ Day(7*h), (x -> x.median[h]).(baseTemp1.forecast[1:end-h]), lw = 2,
      label = "ARMA-l")
plot!(baseTemp2.date[1:end-h] .+ Day(7*h), (x -> x.median[h]).(baseTemp2.forecast[1:end-h]), lw = 2,
      label = "OLS-l")
#
xlims!((minimum(datesTemp), maximum(datesTemp) + Day(7*5)))
ylims!((0, 2.5*10^6))

datesTemp = sort(unique(cfc.date[cfc.model .== 46]))
indTemp = findall((covidTruth.location .== 0) .& ((d -> d in datesTemp).(covidTruth.date)))
plot(covidTruth.date[covidTruth.location .== 0], log.(covidTruth.value[covidTruth.location .== 0] .+ 1),
     label = "Truth", xlabel = "Date", ylabel = "Cases (log)", lw = 3)
cfcTemp = cfc[(cfc.model .== 46) .& (cfc.location .== 0), :]
plot!(cfcTemp.date .+ Day(7), (x -> log(x.median[1] .+ 1)).(cfcTemp.forecast), lw = 2,
      label = "h = 1")
plot!(baseTemp1.date .+ Day(7), (x -> log(x.median[1] + 1)).(baseTemp1.forecast), lw = 2, label = "ARMA-l")
plot!(baseTemp2.date .+ Day(7), (x -> log(x.median[1] + 1)).(baseTemp2.forecast), lw = 2, label = "OLS-l")
xlims!((minimum(datesTemp), maximum(datesTemp) + Day(7*5)))
ylims!((11, 16.5))


### Where are models good at forecasting?
# Select only one of the well performing models and vary the baseline

fluEns = deepcopy(ifc[ifc.model .== 11, :])
relWISloc = zeros(52, 21, 2)

for l = 1:52
    temp = deepcopy(fluEns.forecast[fluEns.location .== locations[l]])
    wEns = (h -> WIS.(temp[(x -> x.horizon[end]).(temp) .>= h], h, logTrafo = false)).(1:4)
    wEnsLog = (h -> WIS.(temp[(x -> x.horizon[end]).(temp) .>= h], h, logTrafo = true)).(1:4)
    
    for b = 1:21
        baseTemp = deepcopy(findA(fluEns, iBase[iBase.model .== b, :])[2])
        temp2 = deepcopy(baseTemp.forecast[baseTemp.location .== locations[l]])
        wBase = (h -> WIS.(temp2[(x -> x.horizon[end]).(temp2) .>= h], h, logTrafo = false)).(1:4)

        relWISloc[l, b, 1] = sum(vcat(wEns...))/sum(vcat(wBase...))

        baseTemp = deepcopy(findA(fluEns, iBase[iBase.model .== b, :])[2])
        temp2 = deepcopy(baseTemp.forecast[baseTemp.location .== locations[l]])
        wBase = (h -> WIS.(temp2[(x -> x.horizon[end]).(temp2) .>= h], h, logTrafo = true)).(1:4)

        relWISloc[l, b, 2] = sum(vcat(wEnsLog...))/sum(vcat(wBase...))


        println([l, b])
    end
end


chEns = deepcopy(cfc[cfc.model .== 6, :])
crelWISloc = zeros(52, 21, 2)

for l = 1:52
    temp = deepcopy(chEns.forecast[chEns.location .== locations[l]])
    wEns = (h -> WIS.(temp[(x -> x.horizon[end]).(temp) .>= h], h, logTrafo = false)).(1:4)
    wEnsLog = (h -> WIS.(temp[(x -> x.horizon[end]).(temp) .>= h], h, logTrafo = true)).(1:4)
    for b = 1:21
        baseTemp = deepcopy(findA(chEns, cBase[cBase.model .== b, :])[2])
        temp2 = deepcopy(baseTemp.forecast[baseTemp.location .== locations[l]])
        wBase = (h -> WIS.(temp2[(x -> x.horizon[end]).(temp2) .>= h], h, logTrafo = false)).(1:4)

        crelWISloc[l, b, 1] = sum(vcat(wEns...))/sum(vcat(wBase...))

        baseTemp = deepcopy(findA(chEns, cBase[cBase.model .== b, :])[2])
        temp2 = deepcopy(baseTemp.forecast[baseTemp.location .== locations[l]])
        wBase = (h -> WIS.(temp2[(x -> x.horizon[end]).(temp2) .>= h], h, logTrafo = true)).(1:4)

        crelWISloc[l, b, 2] = sum(vcat(wEnsLog...))/sum(vcat(wBase...))
        println([l, b])
    end
end

oi = sortperm(iDivs2)
oc = sortperm(cDivs2)

# Sorting should be done by location stratified CvM
iDivsLoc = zeros(21)
cDivsLoc = zeros(21)

for b = 1:21
    iDivsLoc[b] = sum((l -> CvMdivergence(iBase.forecast[(iBase.model .== b) .& (iBase.location .== l)], [1, 2, 3, 4])).(locations))
    cDivsLoc[b] = sum((l -> CvMdivergence(cBase.forecast[(cBase.model .== b) .& (cBase.location .== l)], [1, 2, 3, 4])).(locations))
end


bar(10.2:-1:1.2, -iDivsLoc[1:10], orientation = :horizontal, legend = :bottomleft,
    bar_width = 0.4, xlabel = "Calibration", label = "Original", ylabel = "Baseline",
    top_margin = 10Plots.mm, dpi = 300)
bar!(9.8:-1:0.8, -[iDivsLoc[11:17]; missing; iDivsLoc[18:19]],
    orientation = :horizontal, bar_width = 0.4, label = "Log")
scatter!(-iDivs[1:10], 10.2:-1:1.2, label = "Global calibration", color = "black")
scatter!(-[iDivs[11:17]; missing; iDivs[18:19]],9.8:-1:0.8, color = "black", label = "")
ylims!((0, 10.5))
yticks!(1:10, modelNamesShort[10:-1:1])
bar!(10.2:-1:1.2, cDivsLoc[1:10], orientation = :horizontal, color = 1,
    bar_width = 0.4, label = "")
bar!(9.8:-1:0.8, [cDivsLoc[11:17]; missing; cDivsLoc[18:19]], color = 2,
    orientation = :horizontal, bar_width = 0.4, label = "")
scatter!(cDivs[1:10], 10.2:-1:1.2, label = "", color = "black")
scatter!([cDivs[11:17]; missing; cDivs[18:19]],9.8:-1:0.8, color = "black", label = "")
vline!([0], label = "", color = "black")
annotate!(-maximum(iDivsLoc)/2, 11, text("Influenza", 14, :center))
annotate!(maximum(cDivsLoc)/2, 11, text("COVID-19", 14, :center))
xticks!(-2000:1000:2000, ["2000", "1000", "0", "1000", "2000"])
# vline!([cDivsLoc[20]], label = "Individual", color = "black", lw = 2, linestyle = :dash)
# vline!([cDivsLoc[21]], label = "Marginal (all)", color = "black", lw = 2, linestyle = :solid)
savefig(path*"/Plots/GlobalVsUniform.png")

oi = sortperm(iDivsLoc)
oc = sortperm(cDivsLoc)

iLocHeat = corkendall(relWISloc[:, :, 1])[oi, oi]
cLocHeat = corkendall(crelWISloc[:, :, 1])[oc, oc]

iLocHeatLog = corkendall(relWISloc[:, :, 2])[oi, oi]
cLocHeatLog = corkendall(crelWISloc[:, :, 2])[oc, oc]


cLower = 0.0
# modelNamesShort = [modelNamesShort; "Indiv"; "Mgn-All"]

mean(iLocHeat .< 0)
mean(iLocHeatLog .< 0)
mean(cLocHeat .< 0)
mean(cLocHeatLog .< 0)

p1 = heatmap(iLocHeat[:, end:-1:1], grid = false, xaxis = false, yaxis = false,
             clims = (cLower, 1), legend = false, left_margin = 5Plots.mm,
             bottom_margin = 5Plots.mm)
for i = 1:21
    temp = modelNamesShort[findfirst(iDivsLoc .== sort(iDivsLoc)[i])]
    annotate!(0, 22 - i, text(temp, 9, :right))
    annotate!(i, 0, text(temp, 9, :right, rotation = 60))
end
title!("WIS")
plot!([0.5, 21.5], [-3.5, -3.5], label = "", arrow = :closed, color = "black")
plot!([-3.5, -3.5], [21.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(11, -4, text("Uniform Calibration", 9))
annotate!(2, -4, text("Good", 9, :left))
annotate!(20, -4, text("Poor", 9, :right))

annotate!(-6, 10, text("Influenza", 14, rotation = 90))
annotate!(-4.5, 11, text("Uniform Calibration", 9, rotation = 90))
annotate!(-4.5, 20, text("Good", 9, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 9, :left, rotation = 90))

p2 = heatmap(iLocHeatLog[:, end:-1:1], grid = false, xaxis = false, yaxis = false,
             clims = (cLower, 1), legend = false, bottom_margin = 5Plots.mm)
for i = 1:21
    temp = modelNamesShort[findfirst(iDivsLoc .== sort(iDivsLoc)[i])]
    annotate!(0, 22 - i, text(temp, 9, :right))
    annotate!(i, 0, text(temp, 9, :right, rotation = 60))
end
title!("lWIS")
plot!([0.5, 21.5], [-3.5, -3.5], label = "", arrow = :closed, color = "black")
plot!([-3.5, -3.5], [21.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(11, -4, text("Uniform Calibration", 9))
annotate!(2, -4, text("Good", 9, :left))
annotate!(20, -4, text("Poor", 9, :right))

annotate!(-4.5, 11, text("Uniform Calibration", 9, rotation = 90))
annotate!(-4.5, 20, text("Good", 9, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 9, :left, rotation = 90))

p3 = heatmap(cLocHeat[:, end:-1:1], grid = false, xaxis = false, yaxis = false,
             clims = (cLower, 1), legend = false, left_margin = 10Plots.mm,
             bottom_margin = 5Plots.mm)
for i = 1:21
    temp = modelNamesShort[findfirst(cDivsLoc .== sort(cDivsLoc)[i])]
    annotate!(0, 22 - i, text(temp, 9, :right))
    annotate!(i, 0, text(temp, 9, :right, rotation = 60))
end
plot!([0.5, 21.5], [-3.5, -3.5], label = "", arrow = :closed, color = "black")
plot!([-3.5, -3.5], [21.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(11, -4, text("Uniform Calibration", 9))
annotate!(2, -4, text("Good", 9, :left))
annotate!(20, -4, text("Poor", 9, :right))

annotate!(-6, 10, text("COVID-19", 14, rotation = 90))
annotate!(-4.5, 11, text("Uniform Calibration", 9, rotation = 90))
annotate!(-4.5, 20, text("Good", 9, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 9, :left, rotation = 90))

p4 = heatmap(cLocHeatLog[:, end:-1:1], grid = false, xaxis = false, yaxis = false,
             clims = (cLower, 1), legend = false, bottom_margin = 5Plots.mm)
for i = 1:21
    temp = modelNamesShort[findfirst(cDivsLoc .== sort(cDivsLoc)[i])]
    annotate!(0, 22 - i, text(temp, 9, :right))
    annotate!(i, 0, text(temp, 9, :right, rotation = 60))
end
plot!([0.5, 21.5], [-3.5, -3.5], label = "", arrow = :closed, color = "black")
plot!([-3.5, -3.5], [21.5, 0.5], label = "", arrow = :closed, color = "black")

annotate!(11, -4, text("Uniform Calibration", 9))
annotate!(2, -4, text("Good", 9, :left))
annotate!(20, -4, text("Poor", 9, :right))

annotate!(-4.5, 11, text("Uniform Calibration", 9, rotation = 90))
annotate!(-4.5, 20, text("Good", 9, :right, rotation = 90))
annotate!(-4.5, 2, text("Poor", 9, :left, rotation = 90))

dummy_data = fill(missing, 2, 2)
p_colorbar = heatmap(dummy_data, clims = (cLower, 1.0), 
                    colorbar = :right,
                    showaxis = false, 
                    grid = false,
                    framestyle = :none,
                    background_color = :transparent,
                    foreground_color = :transparent,
                    ytickfontsize = 10)
#
pEmpty = plot([1, 2], legend = false, axis = false, grid = false, color = :transparent)

l1 = @layout [a{0.1h}; b{0.8h}; c{0.1h}]
pBar = plot(pEmpty, p_colorbar, pEmpty, layout = l1, size = (100, 900))

# Use a custom layout
l = @layout [[a b; c d] e{0.05w}]
plot(p1, p2, p3, p4, pBar, layout = l, size = (1000, 900), dpi = 300)
savefig(path*"/Plots/LocHeatmap.png")

mean(CI .< 0.9)
mean(CIl .< 0.9)
mean(CC .< 0.9)
mean(CCl .< 0.9)