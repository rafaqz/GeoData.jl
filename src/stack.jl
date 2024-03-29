# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
    AbstractGeoStack

Abstract supertype for objects that hold multiple [`AbstractGeoArray`](@ref)
that share spatial dimensions.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractGeoArray`](@ref), string paths that will load [`AbstractGeoArray`](@ref),
or a single path that points to as a file itself containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`AbstractGeoStack` can hold layers that share some or all of their dimensions.
They cannot have the same dimension with t different length or spatial extent as
another layer.

`getindex` on a `AbstractGeoStack` generally returns a memory backed standard
[`GeoArray`](@ref). `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on a `AbstractGeoStack` with a key returns another stack with
getindex applied to all the arrays in the stack.
"""
abstract type AbstractGeoStack{L} <: AbstractDimStack{L} end

layermissingval(stack::AbstractGeoStack) = stack.layermissingval
filename(stack::AbstractGeoStack) = filename(parent(stack))
missingval(s::AbstractGeoStack, key::Symbol) = _singlemissingval(layermissingval(s), key)

isdisk(A::AbstractGeoStack) = isdisk(first(A))

"""
    subset(s::AbstractGeoStack, keys)

Subset a stack to hold only the layers in `keys`, where `keys` is a `Tuple`
or `Array` of `String` or `Symbol`, or a `Tuple` or `Array` of `Int`
"""
subset(s::AbstractGeoStack, keys) = subset(s, Tuple(keys))
function subset(s::AbstractGeoStack, keys::NTuple{<:Any,<:Key})
    GeoStack(map(k -> s[k], Tuple(keys)))
end
function subset(s::AbstractGeoStack, I::NTuple{<:Any,Int})
    subset(s, map(i -> keys(s)[i], I))
end

_singlemissingval(mvs::NamedTuple, key) = mvs[key]
_singlemissingval(mv, key) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractGeoStack) = invoke(DD.DimTable, Tuple{AbstractDimStack}, read(stack))

function DD.layers(s::AbstractGeoStack{<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(
    s::AbstractGeoStack, data, dims=dims(s), refdims=refdims(s), 
    layerdims=DD.layerdims(s), metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s), 
)
    DD.basetypeof(s)(data, dims, refdims, layerdims, metadata, layermetadata, layermissingval)
end
function DD.rebuild(s::AbstractGeoStack;
    data=parent(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s),
)
    DD.basetypeof(s)(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval
    )
end

function DD.rebuild_from_arrays(
    s::AbstractGeoStack, das::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractDimArray}}}; 
    refdims=DD.refdims(s), 
    metadata=DD.metadata(s), 
    data=map(parent, das), 
    dims=DD.combinedims(das...), 
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    layermissingval=map(missingval, das),
)
    rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, layermissingval)
end

# Base methods #################################################################

Base.names(s::AbstractGeoStack) = keys(s)
Base.copy(stack::AbstractGeoStack) = map(copy, stack)

#### Stack getindex ####
# Different to DimensionalData as we construct a GeoArray
Base.getindex(s::AbstractGeoStack, key::AbstractString) = s[Symbol(key)]
function Base.getindex(s::AbstractGeoStack, key::Symbol)
    data_ = parent(s)[key]
    dims_ = dims(s, DD.layerdims(s, key))
    metadata = DD.layermetadata(s, key)
    GeoArray(data_, dims_, refdims(s), key, metadata, missingval(s, key))
end
# Key + Index
@propagate_inbounds @inline function Base.getindex(s::AbstractGeoStack, key::Symbol, i1, I...)
    A = s[key][i1, I...]
end


# Concrete AbstrackGeoStack implementation #################################################

"""
    GeoStack <: AbstrackGeoStack

    GeoStack(data...; name, kw...)
    GeoStack(data::Union{Vector,Tuple}; name, kw...)
    GeoStack(data::NamedTuple; kw...))
    GeoStack(s::AbstractGeoStack; kw...)
    GeoStack(s::AbstractGeoArray; layersfrom=Band, kw...)

Load a file path or a `NamedTuple` of paths as a `GeoStack`, or convert arguments, a 
`Vector` or `NamedTuple` of `GeoArray` to `GeoStack`.

# Arguments

- `data`: A `NamedTuple` of [`GeoArray`](@ref), or a `Vector`, `Tuple` or splatted arguments
    of [`GeoArray`](@ref). The latter options must pass a `name` keyword argument.

# Keywords

- `name`: Used as stack layer names when a `Tuple`, `Vector` or splat of `GeoArray` is passed in.
- `metadata`: A `Dict` or `DimensionalData.Metadata` object.
- `refdims`: `Tuple` of `Dimension` that the stack was sliced from.
- `layersfrom`: `Dimension` to source stack layers from if the file is not 
    already multi-layered. This will often be `Band`, which is the default.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = GeoStack(files; mappedcrs=EPSG(4326))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
struct GeoStack{L<:Union{FileStack,NamedTuple},D<:Tuple,R<:Tuple,LD<:NamedTuple,M,LM,LMV} <: AbstractGeoStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    layermissingval::LMV
end
# Multi-file stack from strings
function GeoStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    name=map(filekey, filenames), keys=name, kw...
)
    GeoStack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
function GeoStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    crs=nothing, mappedcrs=nothing, source=nothing, kw...
) where K
    layers = map(keys(filenames), values(filenames)) do key, fn
        source = source isa Nothing ? _sourcetype(fn) : source
        crs = defaultcrs(source, crs)
        mappedcrs = defaultmappedcrs(source, mappedcrs)
        _open(fn; key) do ds
            data = FileArray(ds, fn; key)
            dims = DD.dims(ds, crs, mappedcrs)
            md = metadata(ds)
            mv = missingval(ds)
            GeoArray(data, dims; name=key, metadata=md, missingval=mv)
        end
    end
    GeoStack(NamedTuple{K}(layers); kw...)
end
# Multi GeoArray stack from Tuple of AbstractArray
function GeoStack(data::Tuple{Vararg{<:AbstractArray}}, dims::DimTuple; name=nothing, keys=name, kw...)
    return GeoStack(NamedTuple{cleankeys(keys)}(data), dims; kw...)
end
# Multi GeoArray stack from NamedTuple of AbstractArray
function GeoStack(data::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}, dims::DimTuple; kw...)
    # TODO: make this more sophisticated an match dimension length to axes?
    layers = map(data) do A
        GeoArray(A, dims[1:ndims(A)])
    end
    return GeoStack(layers; kw...)
end
# Multi GeoArray stack from AbstractDimArray splat
GeoStack(layers::AbstractDimArray...; kw...) = GeoStack(layers; kw...)
# Multi GeoArray stack from tuple with `keys` keyword
function GeoStack(layers::Tuple{Vararg{<:AbstractGeoArray}}; 
    name=map(name, layers), keys=name, kw...
)
    GeoStack(NamedTuple{cleankeys(keys)}(layers); kw...)
end
# Multi GeoArray stack from NamedTuple
function GeoStack(layers::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}};
    resize=nothing, refdims=(), metadata=NoMetadata(), kw...
)
    # resize if not matching sizes - resize can be `crop`, `resample` or `extend`
    layers = resize isa Nothing ? layers : resize(layers)
    # DD.comparedims(layers...)
    dims=DD.combinedims(layers...)
    data=map(parent, layers)
    layerdims=map(DD.basedims, layers)
    layermetadata=map(DD.metadata, layers)
    layermissingval=map(missingval, layers)
    return GeoStack(
        data, dims, refdims, layerdims, metadata,
        layermetadata, layermissingval
    )
end
# Single-file stack from a string
function GeoStack(filename::AbstractString;
    dims=nothing, refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing,
    layerdims=nothing, layermetadata=nothing, layermissingval=nothing,
    source=_sourcetype(filename), name=nothing, keys=name, layersfrom=Band,
    resize=nothing,
)
    st = if haslayers(_sourcetype(filename))
        crs = defaultcrs(source, crs)
        mappedcrs = defaultmappedcrs(source, mappedcrs)
        data, field_kw = _open(filename) do ds
            dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
            refdims = refdims == () || refdims isa Nothing ? () : refdims
            layerdims = layerdims isa Nothing ? DD.layerdims(ds) : layerdims
            metadata = metadata isa Nothing ? DD.metadata(ds) : metadata
            layermetadata = layermetadata isa Nothing ? DD.layermetadata(ds) : layermetadata
            layermissingval = layermissingval isa Nothing ? GeoData.layermissingval(ds) : layermissingval
            data = FileStack{source}(ds, filename; keys)
            data, (; dims, refdims, layerdims, metadata, layermetadata, layermissingval)
        end
        GeoStack(data; field_kw...)
    else
        # Band dims acts as layers
        GeoStack(GeoArray(filename); layersfrom)
    end

    # Maybe split the stack into separate arrays to remove extra dims.
    if !(keys isa Nothing)
        return map(identity, st)
    else
        return st
    end
end
function GeoStack(A::GeoArray; 
    layersfrom=Band, name=nothing, keys=name, metadata=metadata(A), refdims=refdims(A), kw...
)
    layersfrom = layersfrom isa Nothing ? Band : layersfrom
    keys = keys isa Nothing ? _layerkeysfromdim(A, layersfrom) : keys
    slices = slice(A, layersfrom)
    layers = NamedTuple{Tuple(map(Symbol, keys))}(Tuple(slices))
    GeoStack(layers; refdims=refdims, metadata=metadata, kw...)
end
# Stack from stack, dims args
GeoStack(st::AbstractGeoStack, dims::DimTuple; kw...) = GeoStack(st; dims, kw...)
# Stack from table, dims args
function GeoStack(table, dims::DimTuple; name=_not_a_dimcol(table, dims), keys=name, kw...)
    # TODO use `name` everywhere, not keys
    if keys isa Symbol
        col = Tables.getcolumn(table, keys)
        layers = NamedTuple{(keys,)}((reshape(col, map(length, dims)),))
    else
        layers = map(keys) do k
            col = Tables.getcolumn(table, k)
            reshape(col, map(length, dims))
        end |> NamedTuple{keys}
    end
    GeoStack(layers, dims; kw...)
end

function _layerkeysfromdim(A, dim)
    map(index(A, dim)) do x
        if x isa Number
            Symbol(string(DD.dim2key(dim), "_", x))
        else
            Symbol(x)
        end
    end
end

# Rebuild from internals
function GeoStack(
    data::Union{FileStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, layermissingval) 
    st = GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval
    )
end
# GeoStack from another stack
function GeoStack(s::AbstractDimStack; name=cleankeys(Base.keys(s)), keys=name,
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s)
)
    st = GeoStack(
        data, DD.dims(s), refdims, layerdims, metadata, layermetadata, layermissingval
    )

    # TODO This is a bit of a hack, it should use `formatdims`. 
    return set(st, dims...)
end

Base.convert(::Type{GeoStack}, src::AbstractDimStack) = GeoStack(src)

GeoArray(stack::GeoStack) = cat(values(stack)...; dims=Band([keys(stack)...]))

defaultcrs(T::Type, crs) = crs
defaultcrs(T::Type, ::Nothing) = defaultcrs(T)
defaultcrs(T::Type) = nothing
defaultmappedcrs(T::Type, crs) = crs
defaultmappedcrs(T::Type, ::Nothing) = defaultmappedcrs(T)
defaultmappedcrs(T::Type) = nothing

# Precompile
precompile(GeoStack, (String,))

@deprecate stack(args...; kw...) GeoStack(args...; kw...)
