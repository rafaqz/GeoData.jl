using HDF5

export SMAPseries, SMAPstack, smapseries

const SMAPMISSING = -9999.0f0
const SMAPGEODATA = "Geophysical_Data"
const SMAPCRS = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
const SMAPSIZE = (3856, 1624)
const SMAPDIMTYPES = (X, Y) 

# Dataset wrapper ###############################################################
# Becauase SMAP is just one of many HDF5 formats,
# we wrap it in SMAPhdf5 and SMAPvar wrappers

struct SMAPhdf5{T}
    ds::T
end

missingval(ds::SMAPhdf5) = SMAPMISSING
layermissingval(ds::SMAPhdf5) = SMAPMISSING
layerkeys(ds::SMAPhdf5) = keys(ds)
layersizes(ds::SMAPhdf5, keys) = map(_ -> SMAPSIZE, keys)
filekey(ds::SMAPhdf5, key::Nothing) = first(keys(ds))

function DD.dims(wrapper::SMAPhdf5)
    dataset = parent(wrapper)
    proj = read(HDF5.attributes(HDF5.root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        extent = HDF5.attributes(HDF5.root(dataset)["Metadata/Extent"])
        lonbounds = read(extent["westBoundLongitude"]), read(extent["eastBoundLongitude"])
        latbounds = read(extent["southBoundLatitude"]), read(extent["northBoundLatitude"])
        latvec = read(HDF5.root(dataset)["cell_lat"])[1, :]
        lonvec = read(HDF5.root(dataset)["cell_lon"])[:, 1]
        lonmode = Mapped(
           order=Ordered(),
           span=Irregular(lonbounds),
           sampling=Intervals(Center()),
           crs=SMAPCRS,
           mappedcrs=EPSG(4326)
        )
        latmode = Mapped(
            order=Ordered(ReverseIndex(), ReverseArray(), ForwardRelation()),
            span=Irregular(latbounds),
            sampling=Intervals(Center()),
            crs=SMAPCRS,
            mappedcrs=EPSG(4326),
        )
        (X(lonvec; mode=lonmode), Y(latvec; mode=latmode))
    else
        error("projection $proj not supported")
    end
end

DD.refdims(wrapper::SMAPhdf5, filename) = (_smap_timedim(_smap_timefromfilename(filename)),)

# TODO actually add metadata to the dict
DD.metadata(wrapper::SMAPhdf5) = Metadata{SMAPfile}(Dict())

function DD.layerdims(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    # All dims are the same
    NamedTuple{keys}(map(_ -> SMAPDIMTYPES, keys))
end

function DD.layermetadata(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    NamedTuple{keys}(map(_ -> DD.metadata(ds), keys))
end

Base.keys(ds::SMAPhdf5) = cleankeys(keys(parent(ds)[SMAPGEODATA]))
Base.parent(wrapper::SMAPhdf5) = wrapper.ds
Base.getindex(wrapper::SMAPhdf5, path) = wrapper.ds[path]

_smappath(key::Key) = SMAPGEODATA * "/" * string(key)

struct SMAPvar{T}
    ds::T
end
Base.parent(wrapper::SMAPvar) = wrapper.ds


# GeoArray ######################################################################

function FileArray(ds::SMAPhdf5, filename::AbstractString; key, kw...)
    FileArray(SMAPvar(HDF5DiskArray(ds[_smappath(key)])), filename; key, kw...)
end
function FileArray(var::SMAPvar, filename::AbstractString; key, kw...)
    T = eltype(parent(var))
    N = length(SMAPSIZE)
    eachchunk = DA.eachchunk(parent(var))
    haschunks = DA.haschunks(parent(var))
    FileArray{SMAPfile,T,N}(filename, SMAPSIZE; key, eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::FileArray{SMAPfile}; kw...)
    _read(var -> f(HDF5DiskArray(var)), SMAPfile, filename(A); key=key(A), kw...)
end
    

# Stack ########################################################################

@deprecate SMAPstack(args...; kw...) GeoStack(args...; source=SMAPfile, kw...)

function FileStack{SMAPfile}(ds::SMAPhdf5, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? layerkeys(ds) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = HDF5DiskArray(ds[_smappath(key)])
        eltype(var), size(var), DA.eachchunk(var), DA.haschunks(var)
    end
    layertypes = NamedTuple{keys}(map(x->x[1], type_size_ec_hc))
    layersizes = NamedTuple{keys}(map(x->x[2], type_size_ec_hc))
    eachchunk = NamedTuple{keys}(map(x->x[3], type_size_ec_hc))
    haschunks = NamedTuple{keys}(map(x->x[4], type_size_ec_hc))
    FileStack{SMAPfile,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

# Series #######################################################################

"""
    SMAPseries(filenames::AbstractString; kw...)
    SMAPseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)

[`GeoSeries`](@ref) loader for SMAP files and whole folders of files,
organised along the time dimension. Returns a [`GeoSeries`](@ref).

# Arguments

- `filenames`: A `String` path to a directory of SMAP files,
    or a vector of `String` paths to specific files.
- `dims`: `Tuple` containing `Ti` dimension for the series.
    Automatically generated form `filenames` unless passed in.

# Keywords

- `kw`: Passed to `GeoSeries`.
"""
function smapseries(dir::AbstractString; kw...)
    smapseries(joinpath.(dir, filter_ext(dir, ".h5")); kw...)
end
function smapseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)
    if dims isa Nothing
        usedpaths = String[]
        timeseries = []
        errors = []
        for filename in filenames
            try
                t = _smap_timefromfilename(filename)
                push!(timeseries, t)
                push!(usedpaths, filename)
            catch e
                push!(errors, e)
            end
        end
        # Use the first files time dim as a template, but join vals into an array of times.
        timedim = _smap_timedim(timeseries)
    else
        usedpaths = filenames
    end
    # Show errors after all files load, or you can't see them:
    if length(errors) > 0
        println("Some errors thrown during file load: ")
        println.(errors)
    end
    # Get the dims once for the whole series
    dims, metadata =_read(SMAPfile, first(filenames)) do ds
        DD.dims(ds), DD.metadata(ds)
    end
    GeoSeries(usedpaths, (timedim,); child=stack, dims, metadata, kw...)
end

@deprecate SMAPseries(args...; kw...) smapseries(args...; kw...)


# Utils ########################################################################

function _read(f, ::Type{SMAPfile}, filepath::AbstractString; key=nothing, kw...)
    if key isa Nothing
        h5open(ds -> cleanreturn(f(SMAPhdf5(ds))), filepath; kw...)
    else
        h5open(ds -> cleanreturn(f(SMAPhdf5(ds)[_smappath(key)])), filepath)
    end
end

function _smap_timefromfilename(filename::String)
    dateformat = DateFormat("yyyymmddTHHMMSS")
    dateregex = r"SMAP_L4_SM_gph_(\d+T\d+)_"
    datematch = match(dateregex, filename)
    if !(datematch === nothing)
        DateTime(datematch.captures[1], dateformat)
    else
        error("Date/time not correctly formatted in path: $filenampathe")
    end
end

_smap_timedim(t::DateTime) = _smap_timedim(t:Hour(3):t)
_smap_timedim(times::AbstractVector) =
    Ti(times, mode=Sampled(Ordered(), Regular(Hour(3)), Intervals(Start())))



# HDF5 DiskArrays ########################################################################################
# Copied from HDF5Utils as it is not up to date with HDF5 0.14/0.15
mutable struct HDF5DiskArray{T,N,C,D,R} <: AbstractDiskArray{T, N}
    ds::HDF5.Dataset
    cs::C
    lo::NTuple{N,Int}
    hi::NTuple{N,Int}
    cache::Array{T,D}
end

Base.size(x::HDF5DiskArray{T, N}) where {T, N} = size(x.ds)::NTuple{N, Int}

DA.haschunks(x::HDF5DiskArray{<:Any,<:Any,Nothing}) = DA.Unchunked()
DA.haschunks(x::HDF5DiskArray) = DA.Chunked()

DA.eachchunk(x::HDF5DiskArray{<:Any, <:Any, <:DA.GridChunks}) = x.cs

DA.readblock!(x::HDF5DiskArray, aout, r::AbstractUnitRange...) = aout .= x.ds[r...]
DA.writeblock!(x::HDF5DiskArray, v, r::AbstractUnitRange...) = x.ds[r...] = v

const _cache_size = Ref(10 * 1024^2)

function set_cache_size(cache_size)
    _cache_size[] = cache_size
end

get_cache_size() = _cache_size[]

get_cache_size(ds::HDF5.Dataset) = _cache_size[] ÷ sizeof(eltype(ds))

function HDF5DiskArray(ds::HDF5.Dataset)
    cs = try
        disable_dag()
        DA.GridChunks(ds, get_chunk(ds))
        enable_dag()
    catch
        nothing
    end
    T, N, C = eltype(ds), ndims(ds), typeof(cs)
    strides = cumprod(collect(size(ds)))
    D = findlast(strides .< get_cache_size(ds))
    D = min(something(D, 1) + 1, N)
    if D == 1
        R = min(get_cache_size(ds), size(ds, 1))
    else
        R = min(ceil(Int, get_cache_size(ds) / strides[D - 1]), size(ds, D))
    end
    lo, hi = ntuple(zero, N), ntuple(zero, N)
    cache = zeros(T, size(ds)[1:(D - 1)]..., 0)
    HDF5DiskArray{T,N,C,D,R}(ds, cs, lo, hi, cache)
end

@generated function _getindex(x::HDF5DiskArray{T, N, C, D, R}, r::Integer...) where {T, N, C, D, R}
    colons = fill(:(:), D - 1)
    rl = [:(r[$d]) for d in 1:(D - 1)]
    rr = [:(r[$d]) for d in (D + 1):N]
    cond = :(r[$D] < x.lo[$D] || r[$D] > x.hi[$D])
    for d in (D + 1):N
        cond = :($cond || r[$d] != x.lo[$d])
    end
    ex = quote
        @inbounds if $cond
            x.lo, x.hi = r, min.(r .+ $R, size(x))
            x.cache = x[$(colons...), x.lo[$D]:x.hi[$D], $(rr...)]
        end
        @inbounds v = x.cache[$(rl...), r[$D] - x.lo[$D] + 1]
        return v
    end
    return ex
end

Base.getindex(x::HDF5DiskArray, r::CartesianIndex) = _getindex(x, Tuple(r)...)
Base.getindex(x::HDF5DiskArray, r::Integer...) = _getindex(x, r...)
Base.getindex(x::HDF5DiskArray, i::Integer) =  getindex(x, CartesianIndices(x)[i])
Base.getindex(x::HDF5DiskArray{T, 1}, i::Integer) where T = _getindex(x, i)

Base._reshape(x::HDF5DiskArray, dims::NTuple{N, Int}) where N = Base.__reshape((x, IndexStyle(x)), dims)

Base.Array(x::HDF5DiskArray) = read(x.ds)

Base.getindex(x::HDF5DiskArray{T, N}, is::Vararg{Union{AbstractVector, Colon}, N}) where {T, N} = getindex(x.ds, is...)

cleanreturn(A::HDF5DiskArray) = Array(A)
