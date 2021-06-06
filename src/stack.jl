# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
    AbstractGeoStack

Abstract supertype for objects that hold multiple [`AbstractGeoArray`](@ref)
that share spatial bounds.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractGeoArray`](@ref), string paths that will load [`AbstractGeoArray`](@ref),
or a single path that points to as a file itself containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`getindex` on a `AbstractGeoStack` generally returns a memory backed standard
[`GeoArray`](@ref). `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on a `AbstractGeoStack` with a key returns another stack with
getindex applied to all the arrays in the stack.
"""
abstract type AbstractGeoStack{L} <: AbstractDimStack{L} end

window(stack::AbstractGeoStack) = stack.window
layermissingval(stack::AbstractGeoStack) = stack.layermissingval
filename(stack::AbstractGeoStack) = filename(data(stack))
missingval(s::AbstractGeoStack, key::Symbol) = _singlemissingval(layermissingval(s), key)

_singlemissingval(mvs::NamedTuple, key) = mvs[key]
_singlemissingval(mv, key) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractGeoStack) = invoke(DD.DimTable, Tuple{AbstractDimStack}, read(stack))

function DD.layers(s::AbstractGeoStack{<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(s::T;
    data=data(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s), window=window(s),
) where T<:AbstractGeoStack
    DD.basetypeof(T)(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval, window
    )
end

# Base methods #################################################################

Base.names(s::AbstractGeoStack) = keys(s)
Base.copy(stack::AbstractGeoStack) = rebuild(stack; data=map(copy, stack))

#### Stack getindex ####
# Different to DimensionalData as we may have a lazy `window` to apply.
# Symbol key
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol)
    data_ = data(s)[key]
    dims_ = dims(s, DD.layerdims(s, key))
    metadata = DD.layermetadata(s, key)
    A = GeoArray(data_, dims_, refdims(s), key, metadata, missingval(s, key))
    window(s) isa Nothing ? A : view(A, window(s)...)
end
# Multilayer indexing
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, i1::Int, I::Int...)
    window(s) isa Nothing ? A[i1, I...] : map(A -> view(A, window(s)...)[i1, I...], data(s))
end
# Key + Index
@propagate_inbounds @inline function Base.getindex(s::AbstractGeoStack, key::Symbol, i1, I...)
    s[key][i1, I...]
end

@propagate_inbounds function Base.view(s::AbstractGeoStack, I...)
    win = if window isa Nothing
        Base.reindex(DD.dims2indices(s, window(s)), DD.dims2indices(s, I))
    else
        I
    end
    rebuild(s; window=win)
end

# Concrete AbstrackGeoStack implementation #################################################

"""
    GeoStack <: AbstrackGeoStack

    GeoStack(data...; keys, kwargs...)
    GeoStack(data::Union{Vector,Tuple}; keys, kwargs...)
    GeoStack(data::NamedTuple; window=nothing, metadata=NoMetadata(), refdims=()))
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

A concrete `AbstractGeoStack` implementation. Holds layers of [`GeoArray`](@ref).

# Arguments

- `data`: A `NamedTuple` of [`GeoArray`](@ref), or a `Vector`, `Tuple` or splatted arguments
    of [`GeoArray`](@ref). The latter options must pass a `keys` keyword argument.

# Keywords

- `keys`: Used as stack keys when a `Tuple` or `Vector` or splat of geoarrays are passed in.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
    contained arrays when they are accessed.
- `refdims`: Reference dimensions from earlier subsetting.
- `metadata`: A `DimensionalData.Metadata` object.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.
"""
struct GeoStack{L,D,R,LD,M,LM,LMV,W} <: AbstractGeoStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    layermissingval::LMV
    window::W
end
# Multi-file stack from strings
function GeoStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    keys=map(filekey, filenames), kw...
)
    GeoStack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
function GeoStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    crs=nothing, mappedcrs=nothing, source=nothing, kw...
) where K
    layers = map(keys(filenames), values(filenames)) do key, fn
        source = source isa Nothing ? _sourcetype(fn) : source
        crs = defaultcrs(source, crs)
        mappecrs = defaultmappedcrs(source, mappedcrs)
        _read(fn; key) do ds
            data = FileArray(ds, fn; key)
            dims = DD.dims(ds, crs, mappedcrs)
            md = metadata(ds)
            mv = missingval(ds)
            GeoArray(data, dims; name=key, metadata=md, missingval=mv)
        end
    end
    GeoStack(NamedTuple{K}(layers); kw...)
end
GeoStack(layers::AbstractDimArray...; kw...) = GeoStack(layers; kw...)
function GeoStack(data::Tuple{Vararg{<:AbstractGeoArray}}; keys=map(name, data), kw...)
    GeoStack(NamedTuple{cleankeys(keys)}(data); kw...)
end
function GeoStack(layers::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}};
    resize=nothing, refdims=(), metadata=NoMetadata(), window=nothing, kw...
)
    # resize if not matching sizes - resize can be `crop` or `extent`
    layers = resize isa Nothing ? layers : resize(layers)
    # DD.comparedims(layers...)
    dims=DD.combinedims(layers...)
    data=map(parent, layers)
    layerdims=map(DD.basedims, layers)
    layermetadata=map(DD.metadata, layers)
    layermissingval=map(missingval, layers)
    GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata,
        layermissingval, window
    )
end
# Single-file stack from a string
function GeoStack(filename::AbstractString;
    dims=nothing, refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing, 
    layerdims=nothing, layermetadata=nothing, layermissingval=nothing,
    source=_sourcetype(filename), keys=nothing, window=nothing
)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    data, field_kw = _read(filename) do ds
        dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
        refdims = refdims == () || refdims isa Nothing ? DD.refdims(ds, filename) : refdims
        layerdims = layerdims isa Nothing ? DD.layerdims(ds) : layerdims
        metadata = metadata isa Nothing ? DD.metadata(ds) : metadata
        layermetadata = layermetadata isa Nothing ? DD.layermetadata(ds) : layermetadata
        layermissingval = layermissingval isa Nothing ? GeoData.layermissingval(ds) : layermissingval
        data = FileStack{source}(ds, filename; keys)
        data, (; dims, refdims, layerdims, metadata, layermetadata, layermissingval)
    end
    GeoStack(data; field_kw..., window)
end
function GeoStack(
    data::Union{FileStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, layermissingval, window=nothing,) GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata,
        layermissingval, window,
    )
end
function GeoStack(s::AbstractDimStack; keys=cleankeys(Base.keys(s)),
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s), window=nothing
)
    GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval, window
    )
end

DD.maybestack(As::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}}) = GeoStack(As)

Base.convert(::Type{GeoStack}, src::AbstractDimStack) = GeoStack(src)

defaultcrs(T::Type, crs) = crs
defaultcrs(T::Type, ::Nothing) = defaultcrs(T)
defaultcrs(T::Type) = nothing
defaultmappedcrs(T::Type, crs) = crs
defaultmappedcrs(T::Type, ::Nothing) = defaultmappedcrs(T)
defaultmappedcrs(T::Type) = nothing
