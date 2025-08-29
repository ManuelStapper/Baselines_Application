using Plots, Shapefile, GeoDataFrames, CSV, DataFrames, LibGEOS, StatsBase, GeoInterface, ColorSchemes

# Set path here
# path = 
shp = Shapefile.Table(path*"Data/Shapefile/StatesZoom.shp")
shpLines = Shapefile.Table(path*"Data/Shapefile/StatesLines.shp")
shpTiles = Shapefile.Table(path*"Data/Shapefile/StatesTiles.shp")
shpHex = Shapefile.Table(path*"Data/Shapefile/StatesHex.shp")

centres = (i -> [mean((z -> z.x).(shp.geometry[i].points)), mean((z -> z.y).(shp.geometry[i].points))]).(1:52)
centre = [mean((x -> x[1]).(centres)), mean((x -> x[2]).(centres))]
dists = (i -> sqrt(sum((centres[i] .- centre).^2))).(1:52)
dists = dists ./ maximum(dists)
cols = ColorSchemes.viridis[ceil.(Int64, dists*256)]


plot(shp.geometry[shp.location .== sort(unique(shp.location))[1]], color = cols[1], grid = false, xaxis = false, yaxis = false)
for i = 2:52
    plot!(shp.geometry[shp.location .== sort(unique(shp.location))[i]], color = cols[i])
end
plot!(shpLines.geometry, color = :transparent)

plot(shpTiles.geometry[1], color = cols[1], aspect_ratio = 1, grid = false, xaxis = false, yaxis = false)
for i = 2:52
    plot!(shpTiles.geometry[i], color = cols[i])
end
annotate!(shpTiles.xCen, shpTiles.yCen, Plots.text.(shpTiles.abbr, color = "white", 12))

plot(shpHex.geometry[1], color = cols[1], aspect_ratio = 1, grid = false, xaxis = false, yaxis = false)
for i = 2:52
    plot!(shpHex.geometry[i], color = cols[i])
end
annotate!(shpHex.xCen, shpHex.yCen, Plots.text.(shpHex.abbr, color = "white", 12))

