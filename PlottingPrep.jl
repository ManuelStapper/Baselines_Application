
#############################################
### File to prepare plotting of US states ###
#############################################

# - Made specifically for the Baselines project
# - Intended to run once and results saved to data folder

using Plots, Shapefile, GeoDataFrames, CSV, DataFrames, LibGEOS, StatsBase, GeoInterface

#############################
### Some helper functions ###
#############################

function selectParts(x::Shapefile.Polygon, indKeep)
    starts = x.parts .+ 1
    stops = [x.parts[2:end]; length(x.points)]
    nParts = length(starts)
    ind = setdiff(1:nParts, indKeep)
    ind = ind[1 .<= ind .<= nParts]

    n = (stops .- starts) .+ 1
    partsOut = Int32[]
    pointsOut = Shapefile.Point[]

    for i = 1:nParts
        if !(i in ind)
            push!(partsOut, length(pointsOut))
            append!(pointsOut, x.points[starts[i]:stops[i]])
        end
    end
    push!(partsOut, length(pointsOut))

    xs = (z -> z.x).(pointsOut)
    ys = (z -> z.y).(pointsOut)
    MBRout = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), minimum(ys))

    Shapefile.Polygon(MBRout, partsOut, pointsOut, x.indexcache)
end

function scaleGeom(x::Shapefile.Polygon, factor::Float64, ref::Shapefile.Rect)
    cenX = (ref.right + ref.left)/2
    cenY = (ref.top + ref.bottom)/2
    xs = (z -> z.x).(x.points)
    ys = (z -> z.y).(x.points)

    xs = (xs .- cenX) .* factor .+ cenX
    ys = (ys .- cenY) .* factor .+ cenY
    
    outPoints = Shapefile.Point.(xs, ys)
    outMBR = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    Shapefile.Polygon(outMBR, x.parts, outPoints, x.indexcache)
end

function scaleGeom(x::Shapefile.Polygon, factor::Float64)
    cenX = (x.MBR.right + x.MBR.left)/2
    cenY = (x.MBR.top + x.MBR.bottom)/2
    xs = (z -> z.x).(x.points)
    ys = (z -> z.y).(x.points)

    xs = (xs .- cenX) .* factor .+ cenX
    ys = (ys .- cenY) .* factor .+ cenY
    
    outPoints = Shapefile.Point.(xs, ys)
    outMBR = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    Shapefile.Polygon(outMBR, x.parts, outPoints, x.indexcache)
end

function shiftGeom(x::Shapefile.Polygon, delta::Vector{Float64})
    xs = (z -> z.x).(x.points)
    ys = (z -> z.y).(x.points)
    xs = xs .+ delta[1]
    ys = ys .+ delta[2]
    outPoints = Shapefile.Point.(xs, ys)
    outMBR = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    Shapefile.Polygon(outMBR, x.parts, outPoints, x.indexcache)
end

function mergeMBR(geoms, ind, returnRect = false)
    mbrs = (x -> x.MBR).(geoms[ind])
    left = mbrs[1].left
    right = mbrs[1].right
    bottom = mbrs[1].bottom
    top = mbrs[1].top
    for i = 2:length(mbrs)
        left = minimum([left, mbrs[i].left])
        right = maximum([right, mbrs[i].right])
        top = maximum([top, mbrs[i].top])
        bottom = minimum([bottom, mbrs[i].bottom])
    end

    xs = [left, right, right, left]
    ys = [top, top, bottom, bottom]
    if returnRect
        return Shapefile.Rect(left, bottom, right, top)
    end
    Shapefile.Polygon(Shapefile.Rect(left, bottom, right, top), [0], Shapefile.Point.(xs, ys))
end

# Translate back to shapefile...
function convertPolygon(x::LibGEOS.Polygon)
    xCoords = GeoInterface.coordinates(x)[1]
    xs = (z -> z[1]).(xCoords)
    ys = (z -> z[2]).(xCoords)
    xPoints = Shapefile.Point.(xs, ys)
    xMBR = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    Shapefile.Polygon(xMBR, [0], xPoints)
end

function convertPolygon(x::LibGEOS.MultiPolygon)
    geoms = getGeometries(x)
    polys = convertPolygon.(geoms)
    mbrTemp = Shapefile.Rect(minimum((z -> z.MBR.left).(polys)),
                       minimum((z -> z.MBR.bottom).(polys)),
                       maximum((z -> z.MBR.right).(polys)),
                       maximum((z -> z.MBR.top).(polys)))
    ptsTemp = (z -> z.points).(polys)
    partsTemp = [0; cumsum(length.(ptsTemp))[1:end-1]]
    ptsTemp = vcat(ptsTemp...)
    Shapefile.Polygon(mbrTemp, partsTemp, ptsTemp)
end

# Reading in data to match locations

cd(path * "/Influenza/FluSightNew/auxiliary-data")
locTable = CSV.read("locations.csv", DataFrame)

# Reading in US states shapefile and restricting to states
# considered in the study
cd(path*"/Shapefile")
shp = Shapefile.Table("us-state-boundaries.shp")
keep = findall(parse.(Int64, shp.state) .<= 56)

# Sort data by the location number
# Creating new shapefile data frame to be saved later
o = sortperm(shp.state[keep])
shp_df = DataFrame(location = shp.state[keep][o], name = shp.name[keep][o], geometry = shp.geometry[keep][o])


# Remove small islands from Alaska and Hawaii
shp_df.geometry[2] = selectParts(shp_df.geometry[2], 1:6)
shp_df.geometry[12] = selectParts(shp_df.geometry[12], 1:5)

# Move and scale Alaska and Hawaii for plotting
shp_df.geometry[2] = shiftGeom(scaleGeom(shp_df.geometry[2], 0.3), [34, -27.0])
shp_df.geometry[12] = shiftGeom(scaleGeom(shp_df.geometry[12], 1.2), [49, 6.0])

# Assing the US main land as additional "state"
# Trick to avoid weird boundaries inbetween: Small buffer around states before merging
usGeom = LibGEOS.buffer(shp_df.geometry[1], 0.01)
for i = setdiff(3:51, 12)
    usGeom = LibGEOS.union(usGeom, LibGEOS.buffer(shp_df.geometry[i], 0.01))
end

usGeom = getGeometries(usGeom)[1]
usCoord = GeoInterface.coordinates(usGeom)[1]

xs = (z -> z[1]).(usCoord)
ys = (z -> z[2]).(usCoord)

usPoints = Shapefile.Point.(xs, ys)
usMBR = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
usGeom = Shapefile.Polygon(usMBR, [0], usPoints)

usGeom = shiftGeom(scaleGeom(usGeom, 0.15), [7, -11.0])
shp_df = [DataFrame(location = "00", name = "US", geometry = usGeom); shp_df]


# Prepare prettier plotting by zoom windows around small states
# "District of Columbia" (10)
# "Rhode Island" (41)
# "Delaware" (9)

# Rough outline:

# plot(shp_df.geometry, grid = false, xaxis = false, yaxis = false,
#      color = "white", legend = false)
temp2 = mergeMBR(shp_df.geometry, [10, 9])
temp2 = scaleGeom(temp2, 1.2)
# plot!(temp2, color = :transparent)
temp2Scl = shiftGeom(scaleGeom(temp2, 4.0), [5, -10.0])
# plot!(temp2Scl, color = :transparent)
# plot!([temp2.MBR.left, temp2Scl.MBR.left], [temp2.MBR.bottom, temp2Scl.MBR.bottom], color = "black")
# plot!([temp2.MBR.right, temp2Scl.MBR.right], [temp2.MBR.top, temp2Scl.MBR.top], color = "black")

temp2 = mergeMBR(shp_df.geometry, [41])
temp2 = scaleGeom(temp2, 1.2)
# plot!(temp2, color = :transparent)
temp2Scl = shiftGeom(scaleGeom(temp2, 4.0), [-7.5, 5.0])
# plot!(temp2Scl, color = :transparent)

# plot!([temp2.MBR.left, temp2Scl.MBR.left], [temp2.MBR.bottom, temp2Scl.MBR.bottom], color = "black")
# plot!([temp2.MBR.right, temp2Scl.MBR.right], [temp2.MBR.top, temp2Scl.MBR.top], color = "black")

shp_df.location = parse.(Int64, shp_df.location)

# Fill the zoom windows with districts:

# Bottom window
temp2 = mergeMBR(shp_df.geometry, [10, 9])
temp2 = scaleGeom(temp2, 1.2)
temp2 = Shapefile.Polygon(temp2.MBR, temp2.parts, [temp2.points; temp2.points[1]])

# Which states overlap with the rectangle?
# Why does it need to run 3 times?!
overlaps = (x -> LibGEOS.intersects(temp2, LibGEOS.buffer(x, 0.0))).(shp_df.geometry)
# Index of shapefile that overlaps
ind = findall(overlaps)
# [9, 10, 22, 32, 40, 48]

overlaps = (i -> LibGEOS.intersection(temp2, LibGEOS.buffer(shp_df.geometry[i], 0.0))).(ind)
overlaps = convertPolygon.(overlaps)
for i = 1:length(overlaps)
    overlaps[i] = shiftGeom(scaleGeom(overlaps[i], 4.0, temp2.MBR), [5, -10.0])
end

# add them to the shapefile
shp_df.name = Vector{String}(shp_df.name)
shp_df = [shp_df; DataFrame(location = shp_df.location[ind], name = shp_df.name[ind] .* "Zoom", geometry = overlaps)]

### Top window
temp2 = mergeMBR(shp_df.geometry, [41])
temp2 = scaleGeom(temp2, 1.2)
temp2 = Shapefile.Polygon(temp2.MBR, temp2.parts, [temp2.points; temp2.points[1]])

# Which states overlap with the rectangle?
overlaps = (x -> LibGEOS.intersects(temp2, LibGEOS.buffer(x, 0.0))).(shp_df.geometry)
ind = findall(overlaps)
# [8, 23, 34, 41]

overlaps = (i -> LibGEOS.intersection(temp2, LibGEOS.buffer(shp_df.geometry[i], 0.0))).(ind)
overlaps = convertPolygon.(overlaps)

# Edit below
for i = 1:length(overlaps)
    overlaps[i] = shiftGeom(scaleGeom(overlaps[i], 4.0, temp2.MBR), [-7.5, 5.0])
end

# add them to the shapefile
shp_df = [shp_df; DataFrame(location = shp_df.location[ind], name = shp_df.name[ind] .* "Zoom", geometry = overlaps)]
# plot(shp_df.geometry, grid = false, xaxis = false, yaxis = false)

# Shapefile.write("StatesZoom.shp", shp_df)

# Create Shapefiles for lines
# Zoom and something around Alaska, Hawaii and US

temp2 = mergeMBR(shp_df.geometry, [10, 9])
temp2 = scaleGeom(temp2, 1.2)
temp2 = Shapefile.Polygon(temp2.MBR, temp2.parts, [temp2.points; temp2.points[1]])
temp2Scl = shiftGeom(scaleGeom(temp2, 4.0), [5, -10.0])

# Top left then clockwise 5 points!

ptsTemp = [temp2.points[1:2]; temp2Scl.points[2:4]; temp2.points[4]; temp2.points[1]]
temp2Hull = Shapefile.Polygon(mergeMBR([temp2, temp2Scl], 1:2, true), [0], ptsTemp)

shpLines = DataFrame(name = ["BottomSmall", "BottomZoom", "BottomHull"],
                     geometry = [temp2, temp2Scl, temp2Hull])
#

# plot(shp_df.geometry, color = :transparent)
# plot!(shpLines.geometry, color = :transparent)

temp2 = mergeMBR(shp_df.geometry, [41])
temp2 = scaleGeom(temp2, 1.2)
temp2 = Shapefile.Polygon(temp2.MBR, temp2.parts, [temp2.points; temp2.points[1]])
temp2Scl = shiftGeom(scaleGeom(temp2, 4.0), [-7.5, 5.0])

# Top left then clockwise 5 points!

# plot(shp_df.geometry, color = :transparent)
# plot!(temp2, color = :transparent)
# plot!(temp2Scl, color = :transparent)

ptsTemp = [temp2Scl.points[1:2]; temp2.points[2:4]; temp2Scl.points[4]; temp2Scl.points[1]]
temp2Hull = Shapefile.Polygon(mergeMBR([temp2Scl, temp2], 1:2, true), [0], ptsTemp)

shpLines = [shpLines; DataFrame(name = ["TopSmall", "TopZoom", "TopHull"],
                     geometry = [temp2, temp2Scl, temp2Hull])]
#
# plot(shp_df.geometry, color = :transparent)
# plot!(shpLines.geometry, color = :transparent)

shp_df.geometry[3] = shiftGeom(shp_df.geometry[3], [0.0, -0.5])
shp_df.geometry[1] = shiftGeom(shp_df.geometry[1], [0.0, -0.2])
shp_df.geometry[1] = shiftGeom(shp_df.geometry[1], [-1, 0.0])

# Now boundary lines for Alaska, Hawaii and US
# plot(shp_df.geometry, color = :transparent, legend = false)
# plot!(shpLines.geometry, color = :transparent)

ptsTemp = Shapefile.Point.([-129, -114, -114, -107, -129, -129], [30.5, 30.5, 28, 23.5, 23.5, 30.5])
mbrTemp = Shapefile.Rect(-129, 23.5, -107, 30.5)
poly1 = Shapefile.Polygon(mbrTemp, [0], ptsTemp)

ptsTemp = Shapefile.Point.([-114, -108.44444, -103, -103, -107, -114, -114], [30.5, 30.5, 27, 23.5, 23.5, 28, 30.5])
mbrTemp = Shapefile.Rect(-114, 23.5, -103, 30.5)
poly2 = Shapefile.Polygon(mbrTemp, [0], ptsTemp)

ptsTemp = Shapefile.Point.([-95.5, -84.5, -84.5, -95.5, -95.5], [28, 28, 23.5, 23.5, 28])
mbrTemp = Shapefile.Rect(-95.5, 23.5, -84.5, 28)
poly3 = Shapefile.Polygon(mbrTemp, [0], ptsTemp)

shpLines = [shpLines; DataFrame(name = ["BoxAlaska", "BoxHawaii", "BoxUS"], geometry = [poly1, poly2, poly3])]

# plot(shp_df.geometry, color = :transparent, grid = false, xaxis = false, yaxis = false)
# plot!(shpLines.geometry, color = :transparent)


# Shapefile.write("StatesZoom.shp", shp_df)
# Shapefile.write("StatesLines.shp", shpLines)

# Additionally use tile plot

temp = ["AK	0	0	0	0	0	0	0	0	VT	NH	ME",
    "WA	ID	MT	ND	MN	WI	IL	MI	NJ	NY	MA	0",
    "OR	NV	WY	SD	IA	IN	OH	PA	MD	CT	RI	0",
    "CA	UT	CO	NE	MO	KY	WV	VA	DC	DE	0	0",
    "0	AZ	NM	KS	AR	TN	NC	SC	0	0	0	0",
    "0	0	TX	OK	LA	MS	AL	GA	0	0	0	0",
    "HI	0	0	0	0	0	0	0	FL	0	0	US"]
#
mat = fill("", 7, 12)
for i = 1:7
    mat[i, :] = split(temp[i], "\t")
end

polygons = Vector{Shapefile.Polygon}(undef, 52)
xCen = zeros(Int64, 52)
yCen = zeros(Int64, 52)
for i = 1:52
    ind = findfirst(mat .== locTable.abbreviation[i])
    ys = 8 - ind[1] .+ [-1, 1, 1, -1, -1] .* 0.45
    xs = ind[2] .+ [1, 1, -1, -1, 1] .* 0.45
    ptsTemp = Shapefile.Point.(xs, ys)
    mbrTemp = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    polygons[i] = Shapefile.Polygon(mbrTemp, [0], ptsTemp)
    xCen[i] = ind[2]
    yCen[i] = 8 - ind[1]
end

tiles = DataFrame(location = shp_df.location[1:52], abbr = locTable.abbreviation[1:52],
                  xCen = xCen, yCen = yCen, geometry = polygons)
#

Shapefile.write("StatesTiles.shp", tiles)


# Also make a hexagon tile map
cd(path*"/Shapefile")

function makeHex(row, col, scl = 1)
    θ = [collect(30:60:360); 30] ./ 360 .* (2*π)
    xs = (cos.(θ) ./ cos(π/6) ./ 2) .* scl
    ys = (sin.(θ) ./ cos(π/6) ./ 2) .* scl
    xNew = xs .+ col .+ ifelse(iseven(row), 1/2, 0)
    yNew = ys .+ (8 - row) .* cos(π/6)
    return xNew, yNew
end

# Alternative 1: A bit more compact
temp = ["0	0	0	0	0	0	0	0	0	0	0	0",
    "AK	0	0	0	0	0	0	0	VT	NH	ME	0",
    "0	WA	MT	ND	MN	WI	0	MI	0	NY	MA	RI",
    "0	ID	WY	SD	IA	IL	IN	OH	PA	NJ	CT	0",
    "0	OR	NV	CO	NE	MO	KY	WV	VA	MD	DE	0",
    "0	CA	UT	NM	KS	AR	TN	NC	SC	DC	0	0",
    "0	0	0	AZ	OK	LA	MS	AL	GA	0	0	0",
    "HI	0	0	0	TX	0	0	FL	0	0	US	0"]
#

# Alternative 2: Closer to true positions
temp = ["AK	0	0	0	0	0	0	0	0	0	0	ME",
    "0	0	0	0	0	0	0	0	0	VT	NH	0",
    "0	WA	MT	ND	MN	WI	0	MI	0	NY	MA	RI",
    "0	ID	WY	SD	IA	IL	IN	OH	PA	NJ	CT	0",
    "0	OR	NV	CO	NE	MO	KY	WV	VA	MD	DE	0",
    "0	CA	UT	NM	KS	AR	TN	NC	SC	DC	0	0",
    "0	0	0	AZ	OK	LA	MS	AL	GA	0	0	0",
    "HI	0	0	0	TX	0	0	FL	0	0	US	0"]
#

mat = fill("", length(temp), 12)
for i = 1:length(temp)
    mat[i, :] = split(temp[i], "\t")
end

abbr = locTable.abbreviation[1:end-1]
loc = locTable.location[1:end-1]
loc[loc .== "US"] .= "00"
loc = parse.(Int64, loc)

polyVec = Shapefile.Polygon[]
xCen = Float64[]
yCen = Float64[]

for i = 1:length(abbr)
    ind = findfirst(mat .== abbr[i])
    push!(xCen, ind[2] .+ ifelse(iseven(ind[1]), 1/2, 0))
    push!(yCen, (size(mat)[1] - ind[1]) .* cos(π/6))
    xs, ys = makeHex(ind[1], ind[2], 0.925)

    brTemp = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    ptsTemp = Shapefile.Point.(xs, ys)
    push!(polyVec, Shapefile.Polygon(mbrTemp, [0], ptsTemp))
end

hexTiles = DataFrame(location = loc, abbr = abbr,
                     xCen = xCen, yCen = yCen, geometry = polyVec)
#

Shapefile.write("StatesHex.shp", hexTiles)