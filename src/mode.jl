"""
    AbstractProjected <: AbstractSampled

Abstract supertype for projected index modes.
"""
abstract type AbstractProjected{O,Sp,Sa} <: AbstractSampled{O,Sp,Sa} end

# For now we just remove CRS on GPU - it often contains strings
Adapt.adapt_structure(to, m::AbstractProjected) = Sampled(order(m), span(m), sampling(m))

"""
    Projected <: AbstractProjected

    Projected(order, span, sampling, crs, mappedcrs)
    Projected(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs, mappedcrs=nothing)

An [`AbstractSampled`]($DDabssampleddocs) `IndexMode` with projections attached.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs)
with the addition of `crs` and `mappedcrs` fields.

If both `crs` and `mappedcrs` fields contain CRS data (in a `GeoFormat` wrapper
from GeoFormatTypes.jl) the selector inputs and plot axes will be converted
from and to the specified `mappedcrs` projection automatically. A common use case
would be to pass `mappedcrs=EPSG(4326)` to the constructor when loading eg. a GDALarray:

```julia
GDALarray(filename; mappedcrs=EPSG(4326))
```

The underlying `crs` will be detected by GDAL.

If `mappedcrs` is not supplied (ie. `mappedcrs=nothing`), the base index will be
shown on plots, and selectors will need to use whatever format it is in.
"""
struct Projected{O<:Order,Sp<:Span,Sa<:Sampling,PC,MC} <: AbstractProjected{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    crs::PC
    mappedcrs::MC
end
function Projected(;
    order=AutoOrder(), span=AutoSpan(), sampling=DD.AutoSampling(), crs, mappedcrs=nothing
)
    Projected(order, span, sampling, crs, mappedcrs)
end

crs(mode::Projected) = mode.crs
crs(mode::IndexMode) = nothing

mappedcrs(mode::Projected) = mode.mappedcrs
mappedcrs(mode::IndexMode) = nothing

function DD.rebuild(
    g::Projected, order=order(g), span=span(g), sampling=sampling(g),
    crs=crs(g), mappedcrs=mappedcrs(g)
)
    Projected(order, span, sampling, crs, mappedcrs)
end

"""
    Mapped <: AbstractProjected

    Mapped(order, span, sampling, crs, mappedcrs)
    Mapped(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs=nothing, mappedcrs)

An [`AbstractSampled`]($DDabssampleddocs) `IndexMode`, where the dimension index has
been mapped to another projection, usually lat/lon or `EPSG(4326)`.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs) with the addition of
`crs` and `mappedcrs` fields.

The mapped dimension index will be used as for [`Sampled`]($DDsampleddocs),
but to save in another format the underlying `projectioncrs` may be used.
"""
struct Mapped{O<:Order,Sp<:Span,Sa<:Sampling,PC,MC} <: AbstractProjected{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    crs::PC
    mappedcrs::MC
end
function Mapped(;
    order=AutoOrder(), span=AutoSpan(), sampling=DD.AutoSampling(), crs=nothing, mappedcrs
)
    Mapped(order, span, sampling, crs, mappedcrs)
end

crs(mode::Mapped, dim) = crs(mode)
crs(mode::Mapped) = mode.crs

mappedcrs(mode::Mapped, dim) = mappedcrs(mode)
mappedcrs(mode::Mapped) = mode.mappedcrs

function DD.rebuild(
    g::Mapped, order=order(g), span=span(g),
    sampling=sampling(g), crs=crs(g), mappedcrs=mappedcrs(g)
)
    Mapped(order, span, sampling, crs, mappedcrs)
end

"""
    convertmode(dstmode::Type{<:IndexMode}, x)

Convert the dimension mode between `Projected` and `Mapped`.
Other dimension modes pass through unchanged.

This is used to e.g. save a netcdf file to GeoTiff.
"""
convertmode(dstmode::Type{<:IndexMode}, A::AbstractDimArray) =
    rebuild(A, data(A), convertmode(dstmode, dims(A)))
convertmode(dstmode::Type{<:IndexMode}, dims::Tuple) =
    map(d -> convertmode(dstmode, d), dims)
convertmode(dstmode::Type{<:IndexMode}, dim::Dimension) =
    convertmode(dstmode, DD.basetypeof(mode(dim)), dim)
# Non-projected IndexMode modess pass through
convertmode(dstmode::Type, srcmode::Type{<:IndexMode}, dim::Dimension) = dim
# AbstractProjected passes through if it's the same as dstmode
convertmode(dstmode::Type{M}, srcmode::Type{M}, dim::Dimension) where M<:AbstractProjected = dim
# Otherwise AbstractProjected needs ArchGDAL
convertmode(dstmode::Type, srcmode::Type{<:AbstractProjected}, dim::Dimension) =
    error("Load ArchGDAL.jl to convert projected dimensions")
# The rest of these methods are in reprojected.jl as they need ArchGDAL
