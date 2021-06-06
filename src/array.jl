
"""
    AbstractGeoArray <: DimensionalData.AbstractDimArray

Abstract supertype for objects that wrap an array (or location of an array) 
and metadata about its contents. It may be memory or hold a `FileArray`, which
holds the filename, and is only opened when required.

`AbstractGeoArray`s inherit from [`AbstractDimArray`]($DDarraydocs)
from DimensionalData.jl. They can be indexed as regular Julia arrays or with
DimensionalData.jl [`Dimension`]($DDdimdocs)s. They will plot as a heatmap in
Plots.jl with correct coordinates and labels, even after slicing with
`getindex` or `view`. `getindex` on a `AbstractGeoArray` will always return
a memory-backed `GeoArray`.
"""
abstract type AbstractGeoArray{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ###########################################################
"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing
missingval(A::AbstractGeoArray) = A.missingval

filename(A::AbstractGeoArray) = filename(data(A))

cleanreturn(A::AbstractGeoArray) = modify(cleanreturn, A)
cleanreturn(x) = x

"""
    crs(x)

Get the projected coordinate reference system of a `Y` or `X` `Dimension`,
or of the `Y`/`X` dims of an `AbstractGeoArray`.

For [`Mapped`](@ref) mode this may be `nothing` as there may be no projected
coordinate reference system at all.
"""
function crs end
function crs(A::AbstractGeoArray)
    if hasdim(A, Y)
        crs(dims(A, Y))
    elseif hasdim(A, X)
        crs(dims(A, X))
    else
        error("No Y or X dimension, crs not available")
    end
end
crs(dim::Dimension) = crs(mode(dim))

"""
    mappedcrs(x)

Get the mapped coordinate reference system for the `Y`/`X` dims of an array.

In [`Projected`](@ref) mode this is used to convert [`Selector`]($DDselectordocs)
values form the mappedcrs defined projection to the underlying projection, and to
show plot axes in the mapped projection.

In `Mapped` mode this is the coordinate reference system of the index values.
"""
function mappedcrs end
function mappedcrs(A::AbstractGeoArray)
    if hasdim(A, Y)
        mappedcrs(dims(A, Y))
    elseif hasdim(A, X)
        mappedcrs(dims(A, X))
    else
        error("No Y or X dimension, mappedcrs not available")
    end
end
mappedcrs(dim::Dimension) = mappedcrs(mode(dim))

for f in (:mappedbounds, :projectedbounds, :mappedindex, :projectedindex)
    @eval ($f)(A::AbstractGeoArray, dims_) = ($f)(dims(A, dims_))
    @eval ($f)(A::AbstractGeoArray) = ($f)(dims(A))
end

# DimensionalData methods

DD.units(A::AbstractGeoArray) = get(metadata(A), :units, nothing)
# Rebuild all types of AbstractGeoArray as GeoArray
function DD.rebuild(
    A::AbstractGeoArray, data, dims::Tuple, refdims, name,
    metadata, missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end
function DD.rebuild(A::AbstractGeoArray;
    data=data(A), dims=dims(A), refdims=refdims(A), name=name(A),
    metadata=metadata(A), missingval=missingval(A)
)
    rebuild(A, data, dims, refdims, name, metadata, missingval)
end

function DD.DimTable(As::Tuple{<:AbstractGeoArray,Vararg{<:AbstractGeoArray}}...)
    DimTable(DimStack(map(read, As...)))
end

# DiskArrays methods

DiskArrays.eachchunk(A::AbstractGeoArray) = DiskArrays.eachchunk(parent(A))
DiskArrays.haschunks(A::AbstractGeoArray) = DiskArrays.haschunks(parent(A))

# Base methods

Base.parent(A::AbstractGeoArray) = data(A)

"""
    open(f, A::AbstractGeoArray; write=false)

`open` is used to open any `AbstractGeoArray` and do multiple operations
on it in a safe way. It's a shorthand for the unexported `OpenGeoArray`
constructor. The `write` keyword opens the file in write mode so that it
can be altered on disk using e.g. a broadcast.

`f` is a method that accepts a single argument - an `OpenGeoArray` object
which is just an `AbstractGeoArray` that holds an open disk - based object.
Often it will be a `do` block:

```julia
ga = geoarray(filepath)
open(ga; write=true) do A
    A[I...] .*= 2 # A is an `OpenGeoArray` wrapping the opened disk-based object.
    # ...  other things you need to do with the open file
end
```

By using a do block to open file we ensure they are always closed again
after we finish working with them.
"""
function Base.open(f::Function, A::AbstractGeoArray; kw...)
    # Open FileArray to expose the actual dataset object, even inside nested wrappers
    select = FileArray
    ignore = Union{Dict,Set,Base.MultiplicativeInverses.SignedMultiplicativeInverse}
    fa = Flatten.flatten(data(A), select, ignore)
    if fa == ()
        f(GeoArray(data(A), dims(A), refdims(A), name(A), metadata(A), missingval(A)))
    else
        open(fa[1]; kw...) do x
            # Rewrap the opened object where the FileArray was
            d = Flatten.reconstruct(data(A), (x,), select, ignore) 
            f(GeoArray(d, dims(A), refdims(A), name(A), metadata(A), missingval(A)))
        end
    end
end
Base.write(A::T) where T <: AbstractGeoArray = write(filename(A), A)

# Concrete implementation ######################################################

"""
    GeoArray <: AbsractGeoArray

    GeoArray(A::AbstractArray{T,N}, dims::Tuple; kw...)
    GeoArray(A::AbstractArray{T,N}; dims, kw...)
    GeoArray(A::AbstractGeoArray; kw...) =

A generic, memory-backed spatial array type. All [`AbstractGeoArray`](@ref) are
converted to `GeoArray` when indexed or otherwise transformed.

# Keywords

- `data`: can replace the data in an `AbstractGeoArray`
- `dims`: `Tuple` of `Dimension`s for the array.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from,
    defaulting to `()`.
- `name`: `Symbol` name for the array.
- `missingval`: Value reprsenting missing values, defaulting to `missing`.
    can be passed it.
- `metadata`: `ArrayMetadata` object for the array, or `NoMetadata()`.
"""
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <: AbstractGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function GeoArray(A::AbstractArray, dims::Tuple;
    refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing
)
    GeoArray(A, DD.formatdims(A, dims), refdims, name, metadata, missingval)
end
function GeoArray(A::AbstractArray;
    dims, refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing
)
    GeoArray(A, DD.formatdims(A, dims), refdims, name, metadata, missingval)
end
function GeoArray(A::AbstractGeoArray;
    data=Array(data(A)), dims=dims(A), refdims=refdims(A),
    name=name(A), metadata=metadata(A), missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end
function GeoArray(filename::AbstractString; key=nothing, kw...)
    isfile(filename) || error("File not found: $filename")
    _read(filename) do ds
        key = filekey(ds, key)
        GeoArray(ds, filename, key; kw...)
    end
end
function GeoArray(ds, filename::AbstractString, key=nothing;
    crs=nothing, mappedcrs=nothing, dims=nothing, refdims=(),
    name=Symbol(key isa Nothing ? "" : string(key)),
    metadata=metadata(ds), missingval=missingval(ds), write=false,
    source=_sourcetype(filename)
)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
    data = FileArray(ds, filename; key, write)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end

@propagate_inbounds function Base.setindex!(A::GeoArray, x, I::DD.StandardIndices...)
    setindex!(data(A), x, I...)
end

filekey(ds, key) = key
