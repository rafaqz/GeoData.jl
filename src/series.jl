"""
    AbstractGeoSeries <: DimensionalData.AbstractDimensionalArray

Abstract supertype for high-level `DimensionalArray` that hold stacks, arrays,
or the paths they can be loaded from. `GeoSeries` are indexed with dimensions
as with a `AbstractGeoArray`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.
As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This allows syntax like:

```julia
series[Time(Near(DateTime(2001, 1))][:temp][Y(Between(70, 150)), X(Between(-20,20))] |> plot`
```

`GeoSeries` is the only concrete implementation. It includes a `chiltype` field
indicating the constructor used then loading stacks or arrays of any type from disk,
and holds a `kwargs` `NamedTuple` that will be splatted into to the keyword arguments
of the `child` constructor. This gives control over the construction of lazy-loaded
files.
"""
abstract type AbstractGeoSeries{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ####################################################

DD.metadata(A::AbstractGeoSeries) = nothing
DD.name(A::AbstractGeoSeries) = NoName()
DD.label(A::AbstractGeoSeries) = ""

"""
    modify(f, series::AbstractGeoSeries)

Apply function `f` to the data of the child object.
If the child is an `AbstractGeoStack` the function will
be passed on to its child `AbstractGeoArray`s.

`f` must return an idenically sized array.

This method triggers a complete rebuild of all objects,
and disk based objects will be transferred to memory.

This is useful for swapping out array backend for an
entire series to `CuArray` from CUDA.jl to copy data to a GPU,
and potentially other types like `DAarray` from Distributed.jl.
"""
DD.modify(f, A::AbstractGeoSeries) = map(child -> modify(f, child), values(A))

"""
    GeoSeries <: AbstractGeoSeries

    GeoSeries(A::AbstractArray{<:AbstractGeoArray}, dims; kw...)
    GeoSeries(A::AbstractArray{<:AbstractGeoStack}, dims; kw...)
    GeoSeries(filenames::AbstractArray{<:AbstractString}, dims; kw...)

Concrete implementation of [`AbstractGeoSeries`](@ref).
Series hold paths to array or stack files, along some dimension(s).

# Keywords

- `dims` known dimensions. These are usually read from the first file in the series
    and are assumed to be _the same for all stacks/arrays in the series_.
- `refdims`: existing reference dimension/s 
- `child`: constructor of child objects - `geoarray` or `stack`
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N}} <: AbstractGeoSeries{T,N,D,A}
    data::A
    dims::D
    refdims::R
end
function GeoSeries(
    data::Array{T}, dims; refdims=(), child=nothing, window=nothing
) where T<:Union{<:AbstractGeoStack,<:AbstractGeoArray}
    ser=  GeoSeries(data, DD.formatdims(data, dims), refdims)
    if window isa Nothing
        ser
    else
        map(x -> view(x, window...), ser)
    end
end
function GeoSeries(
    data::Array{T}, dims; refdims=(), child=geoarray, kw...
) where T<:Union{<:AbstractString}
    source = _sourcetype(first(data))
    # Load the firs child
    child1 = child(first(data); source, refdims, kw...)
    if child == geoarray
        # We assume all dims, metadata and missingvals are the same
        childdims = DD.dims(child1)
        metadata = DD.metadata(child1)
        missingval = GeoData.missingval(child1)
        data = map(data) do x
            child(x; 
                dims=childdims, source, metadata, missingval, name=name(child1), kw...
            )
        end
    else
        # We assume all dims, metadata and missingvals are the same
        childdims = DD.dims(child1)
        metadata = DD.metadata(child1)
        layerdims = DD.layerdims(child1)
        layermetadata = DD.layermetadata(child1)
        layermissingval = GeoData.layermissingval(child1)
        data = map(data) do x
            child(x; 
                dims=childdims, source, metadata, layerdims, layermetadata, 
                layermissingval, keys=keys(child1), kw...
            )
        end
    end
    GeoSeries(data, DD.formatdims(data, dims), refdims)
end

@inline function DD.rebuild(
    A::GeoSeries, data, dims::Tuple, refdims=(), name=nothing, metadata=nothing,
)
    GeoSeries(data, dims, refdims)
end
@inline function DD.rebuild(
    A::GeoSeries; 
    data=data(A), dims=dims(A), refdims=refdims(A), name=nothing, metadata=nothing,
)
    GeoSeries(data, dims, refdims)
end
