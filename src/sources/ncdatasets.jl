
export NCDstack, NCDarray

const NCD = NCDatasets
const NCD_STACK_METADATA_KEY = :_ncd_stack_metadata_

const UNNAMEDNCDfile_KEY = "unnamed"

const NCD_FILL_TYPES = (Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String)

# CF standards don't enforce dimension names.
# But these are common, and should take care of most dims.
const NCD_DIMMAP = Dict(
    "lat" => Y,
    "latitude" => Y,
    "lon" => X,
    "long" => X,
    "longitude" => X,
    "time" => Ti,
    "lev" => Z,
    "mlev" => Z,
    "level" => Z,
    "vertical" => Z,
    "x" => X,
    "y" => Y,
    "z" => Z,
    "band" => Band,
)

cansavestack(::Type{NCDfile}) = true
defaultcrs(::Type{NCDfile}) = EPSG(4326) 
defaultmappedcrs(::Type{NCDfile}) = EPSG(4326) 

# GeoArray ########################################################################

@deprecate NCDarray(args...; kw...) GeoArray(args...; source=NCDfile, kw...)

function GeoArray(ds::NCD.NCDataset, filename::AbstractString, key=nothing; kw...)
    key = _firstkey(ds, key)
    GeoArray(ds[key], filename, key; kw...)
end

function _firstkey(ds::NCD.NCDataset, key=nothing)
    key = (key isa Nothing) ? first(layerkeys(ds)) : key |> Symbol
end

function FileArray(var::NCD.CFVariable, filename::AbstractString; kw...)
    da = GeoDiskArray(var)
    size_ = size(da)
    eachchunk = DA.eachchunk(da) 
    haschunks = DA.haschunks(da)
    T = eltype(var)
    N = length(size_)
    FileArray{NCDfile,T,N}(filename, size_; eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::FileArray{NCDfile}; write=A.write, kw...)
    _read(NCDfile, filename(A); key=key(A), write, kw...) do var
        f(GeoDiskArray(var, DA.eachchunk(A), DA.haschunks(A)))
    end
end

"""
    Base.write(filename::AbstractString, ::Type{NCDarray}, s::AbstractGeoArray)

Write an NCDarray to a NetCDF file using NCDatasets.jl

Returns `filename`.
"""
function Base.write(filename::AbstractString, ::Type{NCDfile}, A::AbstractGeoArray)
    open(A) do o
        md = metadata(o)
        key = NCD_STACK_METADATA_KEY
        attribvec = if haskey(md, key) && md[key] isa Metadata{NCDfile} 
            [_stringdict(md[key])...]
        else
            []
        end
        ds = NCD.Dataset(filename, "c"; attrib=attribvec)
        try
            _ncdwritevar!(ds, o)
        finally
            close(ds)
        end
    end
    return filename
end

# Stack ########################################################################

@deprecate NCDstack(args...; kw...) GeoStack(args...; source=NCDfile, kw...)

"""
    Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently `Metadata` is not handled for dimensions, and `Metadata`
from other [`AbstractGeoArray`](@ref) @types is ignored.
"""
function Base.write(filename::AbstractString, ::Type{NCDfile}, s::AbstractGeoStack)
    ds = NCD.Dataset(filename, "c"; attrib=_stringdict(metadata(s)))
    try map(key -> _ncdwritevar!(ds, s[key]), keys(s))
    finally
        close(ds)
    end
end 

# DimensionalData methods for NCDatasets types ###############################

function DD.dims(ds::NCD.Dataset, crs=nothing, mappedcrs=nothing)
    map(_dimkeys(ds)) do key
        _ncddim(ds, key, crs, mappedcrs)
    end |> Tuple
end
function DD.dims(ds::NCD.Dataset, key::Key, crs=nothing, mappedcrs=nothing)
    DD.dims(ds[key], crs, mappedcrs)
end
function DD.dims(var::NCD.CFVariable, crs=nothing, mappedcrs=nothing)
    names = NCD.dimnames(var)
    map(names) do name
        _ncddim(var.var.ds, name, crs, mappedcrs)
    end |> Tuple
end

DD.refdims(ds::NCD.Dataset, filename) = ()

DD.metadata(ds::NCD.Dataset) = Metadata{NCDfile}(DD.metadatadict(ds.attrib))
DD.metadata(ds::NCD.Dataset, key::Key) = metadata(ds[string(key)])
DD.metadata(var::NCD.CFVariable) = Metadata{NCDfile}(DD.metadatadict(var.attrib))
DD.metadata(var::NCD.CFVariable, stackmetadata::AbstractMetadata) = begin
    md = Metadata{NCDfile}(DD.metadatadict(var.attrib))
    md[NCD_STACK_METADATA_KEY] = stackmetadata
    md
end

@inline function DD.layerdims(var::NCD.CFVariable)
    map(NCD.dimnames(var)) do dimname
        _ncddimtype(dimname)()
    end
end
@inline function DD.layerdims(ds::NCD.Dataset)
    keys = Tuple(layerkeys(ds))
    NamedTuple{map(Symbol, keys)}(map(k -> DD.layerdims(ds[string(k)]), keys))
end

@inline function DD.layermetadata(ds::NCD.Dataset)
    keys = Tuple(layerkeys(ds))
    NamedTuple{map(Symbol, keys)}(map(k -> DD.metadata(ds[string(k)]), keys))
end

missingval(var::NCD.CFVariable) = missing

layermissingval(ds::NCD.Dataset) = missing

function layerkeys(ds::NCD.Dataset)
    dimkeys = _dimkeys(ds)
    toremove = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = [ds[k].attrib["bounds"] for k in dimkeys if haskey(ds[k].attrib, "bounds")]
        union(dimkeys, boundskeys)
    else
        dimkeys
    end
    setdiff(keys(ds), toremove)
end

layertypes(ds::NCD.NCDataset, keys) = map(k -> eltype(ds[k]), keys)
layersizes(ds::NCD.NCDataset, keys) = map(k -> size(ds[k]), keys)

function FileStack{NCDfile}(ds::NCD.Dataset, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? layerkeys(ds) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = GeoDiskArray(ds[string(key)])
        eltype(var), size(var), DA.eachchunk(var), DA.haschunks(var)
    end
    layertypes = NamedTuple{keys}(map(x->x[1], type_size_ec_hc))
    layersizes = NamedTuple{keys}(map(x->x[2], type_size_ec_hc))
    eachchunk = NamedTuple{keys}(map(x->x[3], type_size_ec_hc))
    haschunks = NamedTuple{keys}(map(x->x[4], type_size_ec_hc))
    FileStack{NCDfile,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

# Utils ########################################################################

function _read(f, ::Type{NCDfile}, filename::AbstractString; key=nothing, write=false)
    mode = write ? "a" : "r"
    if key isa Nothing
        NCD.Dataset(cleanreturn ∘ f, filename, mode)
    else
        NCD.Dataset(filename, mode) do ds
            cleanreturn(f(ds[_firstkey(ds, key)]))
        end
    end
end

cleanreturn(A::NCD.CFVariable) = Array(A)

function _ncddim(ds, dimname::Key, crs=nothing, mappedcrs=nothign)
    if haskey(ds, dimname)
        dvar = ds[dimname]
        dimtype = _ncddimtype(dimname)
        index = dvar[:]
        meta = Metadata{NCDfile}(DD.metadatadict(dvar.attrib))
        mode = _ncdmode(ds, dimname, index, dimtype, crs, mappedcrs, meta)
        dimtype(index, mode, meta)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoIndex`
        len = _ncfinddimlen(ds, dimname)
        len === nothing && _unuseddimerror()
        Dim{Symbol(dimname)}(Base.OneTo(len), NoIndex(), NoMetadata())
    end
end

function _ncfinddimlen(ds, dimname) 
    for key in keys(ds)
        var = ds[key]
        dimnames = NCD.dimnames(var)
        if dimname in dimnames 
            return size(var)[findfirst(==(dimname), dimnames)]
        end
    end
    return nothing
end

# Find the matching dimension constructor. If its an unknown name 
# use the generic Dim with the dim name as type parameter
_ncddimtype(dimname) = haskey(NCD_DIMMAP, dimname) ? NCD_DIMMAP[dimname] : DD.basetypeof(DD.key2dim(Symbol(dimname)))

_ncfilenamekey(filenames) = cleankeys(layerkeys(fn) for fn in filenames)

function _ncdmode(
    ds, dimname, index::AbstractArray{<:Union{Number,Dates.AbstractTime}}, 
    dimtype, crs, mappedcrs, metadata
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
        # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    # Unless its a time dimension.
    order = _ncdorder(index)
    span, sampling = if haskey(ds[dimname].attrib, "bounds")
        boundskey = ds[dimname].attrib["bounds"]
        boundsmatrix = Array(ds[boundskey])
        Explicit(boundsmatrix), Intervals(Center())
    elseif eltype(index) <: Dates.AbstractTime
        _ncdperiod(index, metadata)
    else
        _ncdspan(index, order), Points()
    end
    if dimtype in (X, Y)
        # If the index is regularly spaced and there is no crs
        # then there is probably just one crs - the mappedcrs
        crs = if crs isa Nothing && span isa Regular
            mappedcrs
        else
            crs
        end
        Mapped(order, span, sampling, crs, mappedcrs)
    else
        Sampled(order, span, sampling)
    end
end
_ncdmode(ds, dimname, index, dimtype, crs, mappedcrs, mode) = Categorical()

function _ncdorder(index)
    index[end] > index[1] ? Ordered(ForwardIndex(), ForwardArray(), ForwardRelation()) :
                            Ordered(ReverseIndex(), ReverseArray(), ForwardRelation())
end

function _ncdspan(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index)))
    step = index[2] - index[1]
    for i in 2:length(index) -1
        # If any step sizes don't match, its Irregular
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if DD.isrev(indexorder(order))
                    index[end] - endhalfcell, index[1] + beginhalfcell
                else
                    index[1] - beginhalfcell, index[end] + endhalfcell
                end
            else
                index[1], index[1]
            end
            return Irregular(bounds)
        end
    end
    # Otherwise regular
    return Regular(step)
end

# delta_t and ave_period are not CF standards, but CDC
function _ncdperiod(index, metadata::Metadata{NCDfile})
    if haskey(metadata, :delta_t)
        period = _parse_period(metadata[:delta_t])
        period isa Nothing || return Regular(period), Points()
    elseif haskey(metadata, :avg_period)
        period = _parse_period(metadata[:avg_period])
        period isa Nothing || return Regular(period), Intervals(Center())
    end
    return sampling = Irregular(), Points()
end

function _parse_period(period_str::String)
    regex = r"(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)"
    mtch = match(regex, period_str)
    if mtch isa Nothing
        nothing
    else
        vals = Tuple(parse.(Int, mtch.captures))
        periods = (Year, Month, Day, Hour, Minute, Second)
        if length(vals) == length(periods)
            compound = sum(map((p, v) -> p(v), periods, vals))
            if length(compound.periods) == 1
                return compound.periods[1]
            else
                return compound
            end
        else
            return nothing 
        end
    end
end

_stringdict(metadata) = attrib = Dict(string(k) => v for (k, v) in metadata)

_dimkeys(ds::NCD.Dataset) = keys(ds.dim)

# Add a var array to a dataset before writing it.
function _ncdwritevar!(ds, A::AbstractGeoArray{T,N}) where {T,N}
    A = reorder(A, ForwardIndex()) |> a -> reorder(a, ForwardRelation())
    # Define required dim vars
    for dim in dims(A)
        dimkey = lowercase(string(name(dim)))
        haskey(ds.dim, dimkey) && continue
        NCD.defDim(ds, dimkey, length(dim))
        mode(dim) isa NoIndex && continue

        # Shift index before conversion to Mapped
        dim = _ncdshiftlocus(dim)
        if dim isa Y || dim isa X
            dim = convertmode(Mapped, dim)
        end
        md = metadata(dim)
        attribvec = md isa Metadata{:NCD} ? [_stringdict(md)...] : []
        if span(dim) isa Explicit
            bounds = val(span(dim))
            boundskey = string(dimkey, "_bnds")
            push!(attribvec, "bounds" => boundskey)
            NCD.defVar(ds, boundskey, bounds, ("bnds", dimkey))
        end
        NCD.defVar(ds, dimkey, Vector(index(dim)), (dimkey,); attrib=attribvec)
    end
    # TODO actually convert the metadata types
    attrib = if metadata isa Metadata{:NCD} 
        _stringdict(metadata(A))
    else
        Dict()
    end
    # Remove stack metdata if it is attached
    pop!(attrib, string(NCD_STACK_METADATA_KEY), nothing)
    # Set _FillValue
    if ismissing(missingval(A))
        eltyp = _notmissingtype(Base.uniontypes(T)...)
        fillval = NCD.fillvalue(eltyp)
        A = replace_missing(A, fillval)
        attrib["_FillValue"] = fillval
    elseif missingval(A) isa T
        attrib["_FillValue"] = missingval(A)
    else
        missingval isa Nothing || @warn "`missingval` $(missingval(A)) is not the same type as your data $T."
    end

    key = if string(name(A)) == ""
        UNNAMEDNCDfile_KEY
    else
        string(name(A))
    end

    dimnames = lowercase.(string.(map(name, dims(A))))
    attribvec = [attrib...] 
    var = NCD.defVar(ds, key, eltype(A), dimnames; attrib=attribvec)
    var[:] = data(A)
end

_notmissingtype(::Type{Missing}, next...) = _notmissingtype(next...)
_notmissingtype(x::Type, next...) = x in NCD_FILL_TYPES ? x : _notmissingtype(next...)
_notmissingtype() = error("Your data is not a type that netcdf can store")

_ncdshiftlocus(dim::Dimension) = _ncdshiftlocus(mode(dim), dim)
_ncdshiftlocus(::IndexMode, dim::Dimension) = dim
function _ncdshiftlocus(mode::AbstractSampled, dim::Dimension)
    if span(mode) isa Regular && sampling(mode) isa Intervals
        # We cant easily shift a DateTime value
        if eltype(dim) isa Dates.AbstractDateTime
            if !(locus(dim) isa Center)
                @warn "To save to netcdf, DateTime values should be the interval Center, rather than the $(nameof(typeof(locus(dim))))"
            end
            dim
        else
            DD.shiftlocus(Center(), dim)
        end
    else
        dim
    end
end

_unuseddimerror(dimname) = error("Netcdf contains unused dimension $dimname")

@inline function _eachchunk(var::NCD.CFVariable) 
    chunkmode, chunkvec = NCDatasets.chunking(var)
    chunksize = chunkmode == :chunked ? Tuple(chunkvec) : size(var)
    DA.GridChunks(var, chunksize)
end

@inline function _haschunks(var::NCD.CFVariable) 
    chunkmode, _ = NCDatasets.chunking(var)
    chunkmode == :chunked ? DA.Chunked() : DA.Unchunked()
end
