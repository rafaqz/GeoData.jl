using GeoData, Plots, Statistics, ArchGDAL, NCDatasets, Dates, DataFrames
using GeoData: stack, Between

geturl(url, filename=splitdir(url)[2]) = begin
    isfile(filename) || download(url, filename)
    filename
end

# Load some layers from NetCDF #############################################
ncurl = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
ncfilename = geturl(ncurl, "tos_O1_2001-2002.nc")
ncstack = stack(ncfilename)

# Load the sea surface temperature layer
A = ncstack[:tos]

# Plot the 1st, 4th, 7th and 10th months 
A[Ti(1:3:12)] |> plot
dims(A, Ti)
savefig("tos_4.png")

# Plot the Australia region
A[Ti(Contains(DateTime360Day(2001, 01, 17))), Lat(Between(-50, 0)), Lon(Between(100, 160))] |> plot
savefig("tos_jan_australia.png")

# Plot the mean
mean(A; dims=Ti) |> plot
savefig("mean_tos.png")

# Save the mean
write("mean_tos.nc", mean(A; dims=Ti))

# Plo a transect
s = read(ncstack)
A[X(Contains(20)), Ti(1)] |> plot
savefig("tos_20deg_lat.png")

# Convert the stack to a dataframe
# TODO broken: 
df = DataFrame(read(ncstack))

# Create a GeoSeries from multiple files.
# This uses the same data three times to avoid downloads, you 
# would really use a series of datasets matching the time dimension dates.
timedim = Ti([DateTime360Day(2001, 01, 1), DateTime360Day(2001, 02, 1), DateTime360Day(2001, 03, 1)])
filenames = [ncfilename, ncfilename, ncfilename]
ncseries = series(filenames, (timedim,); child=stack)


# Get a single array from the series
A = ncseries[Near(DateTime360Day(2001, 01, 1))][:tos][Lon(Between(50, 200))]

# Plot single slices
view(A, Ti(4)) |> plot

# Or reduce over some dimension
mean(A; dims=Ti) |> plot
std(A; dims=Ti) |> plot

# Replace missing values and to reduce over a dimension
minimum(replace_missing(A, NaN); dims=Ti) |> plot
maximum(replace_missing(A, NaN); dims=Ti) |> plot
reduce(+, A; dims=Ti) |> plot

# Plot the mean sea surface temperature for australia in the second half of 2002 
ncstack = NCDstack(ncfilename);
t = Ti(Between(DateTime360Day(2002, 07, 1), DateTime360Day(2002, 012, 30)))
ncstack[:tos][t, Lat(Between(-45, 0.5)), Lon(Between(110, 160))] |> x->mean(x; dims=Ti) |> plot

# Permute the dimensions in the underlying data
# It stil plots the right way up
permutedims(A, (Lat, Lon, Ti))[Ti(1:3:12)] |> plot

# Line plots have (kind of) useful labels
A = ncseries[1][:tos]
A[Lon(170), Ti(10)] |> plot
A[Lat(1:80), Lon(170), Ti(10)] |> plot


# Load a tif with GDAL ######################################################
gdal_url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
tif_filename = geturl(gdal_url)

array = geoarray(tif_filename; name=:cea)
array[Band(1)] |> plot

# Make a stack (just using the same file three times for demonstration)
filepaths = (one=tif_filename, two=tif_filename, three=tif_filename)
diskstack = stack(filepaths)
# Plot a section of the tif file
diskstack[:two][Band(1), Lat(1:100), Lon(1:150)] |> plot

# Convert the data and dims to a DataFrame
@time df = DataFrame(array)

# See the test suit for more examples...
