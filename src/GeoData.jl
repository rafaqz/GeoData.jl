module GeoData

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end GeoData

using Adapt,
      ConstructionBase,
      Dates,
      DiskArrays,
      Missings,
      Mmap,
      ProgressMeter,
      RecipesBase,
      Reexport

import Flatten,
       Setfield,
       HDF5,
       NCDatasets,
       ArchGDAL

@reexport using DimensionalData, GeoFormatTypes, RasterDataSources

const DD = DimensionalData
const DA = DiskArrays

using Base: tail, @propagate_inbounds

using DimensionalData: StandardIndices

using Setfield: @set, @set!


export AbstractGeoArray, GeoArray

export AbstractGeoStack, GeoStack

export AbstractGeoSeries, GeoSeries

export Projected, Mapped

export Band, Lat, Lon, Vert, GeoXDim, GeoYDim, GeoZDim

export missingval, boolmask, missingmask, replace_missing,
       aggregate, aggregate!, disaggregate, disaggregate!,
       crop, extend, slice

export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds

export geoarray, stack, series

const Lon = X
const Lat = Y
const Vert = Z
const GeoXDim = XDim
const GeoYDim = YDim
const GeoZDim = ZDim

struct NCDfile end
struct GRDfile end
struct GDALfile end
struct SMAPfile end

# DimensionalData documentation urls
const DDdocs = "https://rafaqz.github.io/DimensionalData.jl/stable/api"
const DDdimdocs = joinpath(DDdocs, "#DimensionalData.Dimension")
const DDarraydocs = joinpath(DDdocs, "#DimensionalData.AbstractDimensionalArray")
const DDabssampleddocs = joinpath(DDdocs, "#DimensionalData.AbstractSampled")
const DDsampleddocs = joinpath(DDdocs, "#DimensionalData.Sampled")
const DDlocusdocs = joinpath(DDdocs, "#DimensionalData.Locus")
const DDselectordocs = joinpath(DDdocs, "#DimensionalData.Selector")
const DDtidocs = joinpath(DDdocs, "#DimensionalData.Ti")

include("mode.jl")
include("dimensions.jl")
include("filearray.jl")
include("array.jl")
include("filestack.jl")
include("stack.jl")
include("series.jl")
include("utils.jl")
include("aggregate.jl")
include("methods.jl")
include("open.jl")
include("read.jl")
include("sources/grd.jl")
include("show.jl")
include("plotrecipes.jl")
include("convenience.jl")

include("sources/smap.jl")
include("sources/ncdatasets.jl")
include("sources/gdal.jl")
include("resample.jl")
include("reproject.jl")
include("sources/rasterdatasources.jl")

end
