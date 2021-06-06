"""
    replace_missing(a::AbstractGeoArray, newmissingval)
    replace_missing(a::AbstractGeoStack, newmissingval)

Replace missing values in the array or stack with a new missing value,
also updating the `missingval` field/s.

A `GeoArray` containing a newly allocated `Array` is always returned,
even when the missing value matches the current value.
"""
function replace_missing(A::AbstractGeoArray, newmissingval=missing)
    A = read(A)
    newdata = if ismissing(missingval(A))
        if ismissing(newmissingval)
            copy(parent(read(A)))
        else
            collect(Missings.replace(parent(A), newmissingval))
        end
    else
        replace(parent(A), missingval(A) => newmissingval)
    end
    rebuild(A; data=newdata, missingval=newmissingval)
end
function replace_missing(stack::AbstractGeoStack, newmissingval=missing)
    rebuild(stack, map(a -> replace_missing(a, newmissingval, values(stack))))
end

"""
    boolmask(A::AbstractArray, [missingval])

Create a mask array of `Bool` values, from any AbstractArray. For `AbstractGeoArray`
the default `missingval` is `missingval(A)`, for all other `AbstractArray`s
it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function boolmask end
boolmask(A::AbstractGeoArray) =
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name=:Bool_mask)
boolmask(A::AbstractArray, missingval::Missing=missing) = (a -> !ismissing(a)).(parent(A))
boolmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> !isnan(a)).(parent(A))
    else
        (a -> !isapprox(a, missingval)).(parent(A))
    end

"""
    missingmask(A::AbstractArray, [missingval])

Create a mask array of `missing` or `true` values, from any AbstractArray.
For `AbstractGeoArray` the default `missingval` is `missingval(A)`,
for all other `AbstractArray`s it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function missingmask end
missingmask(A::AbstractGeoArray) =
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name=:Missing_mask)
missingmask(A::AbstractArray, missingval::Missing=missing) =
    (a -> ismissing(a) ? missing : true).(parent(A))
missingmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> isnan(a) ? missing : true).(parent(A))
    else
        (a -> isapprox(a, missingval) ? missing : true).(parent(A))
    end

"""
    crop(A::AbstractGeoArray...)

Crop multiple [`AbstractGeoArray`](@ref) to match the size 
of the smallest one for any dimensions that are shared.
"""
crop(layers::NamedTuple{K}) where K = NamedTuple{K}(crop(layers...))
function crop(layers::AbstractGeoArray...)
    dims = DD.combinedims(layers...; check=false)
    alldims = map(DD.dims, layers)
    smallestdims = map(dims) do d
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        reduce(matchingdims) do a, b 
            _choose(_shortest, a, b)
        end |> first
    end
    selectors = map(smallestdims) do sm
        DD.basetypeof(sm)(Between(DD.bounds(sm)))
    end
    map(l -> view(l, selectors...), layers)
end

_choose(f, ::Tuple{}, ::Tuple{}) = ()
_choose(f, ::Tuple{}, (b,)::Tuple) = (b,)
_choose(f, (a,)::Tuple, ::Tuple{}) = (a,) 
_choose(f, (a,)::Tuple, (b,)::Tuple) = (f(a, b) ? a : b,)

_shortest(a, b) = length(a) <= length(b)
_longest(a, b) = length(a) >= length(b)

"""
    extend(A::AbstractGeoArray...)
    extend(A::AbstractGeoArray, dims::Tuple)

Extend multiple [`AbstractGeoArray`](@ref) to match
the size of the largest one. A single `AbstractGeoArray` 
can be extended bu passing the new `dims` tuple as the second
argumen.
"""
extend(layers::NamedTuple{K}) where K = NamedTuple{K}(extend(layers...))
function extend(layers::AbstractGeoArray...)
    dims = DD.combinedims(layers...; check=false) 
    alldims = map(DD.dims, layers)
    largestdims = map(dims) do d
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        reduce(matchingdims) do a, b 
            _choose(_longest, a, b)
        end |> first
    end
    map(l -> extend(l, largestdims), layers)
end
function extend(A::AbstractGeoArray, newdims::Tuple)
    size = map(length, newdims)
    elt = eltype(A)
    newdata = similar(parent(A), elt, size)
    newdata .= missingval(A)
    newA = rebuild(A; data=newdata, dims=newdims)
    ranges = map(dims(A), newdims) do d, nd
        (DD.sel2indices(nd, Near(first(d)))):(DD.sel2indices(nd, Near(last(d))))
    end
    copyto!(parent(newA), CartesianIndices((ranges...,)), 
            parent(read(A)), CartesianIndices(A)) 
    newA
end

"""
    slice(A::Union{AbstractGeoArray,AbstractGeoStack,AbstracGeoSeries}, dims)

Slice an object allong some dimension(s), lazily using `view`. For a single `GeoArray` 
or `GeoStack` this will return a `GeoSeries` of `GeoArray` or `GeoStack` that are slices 
along the specified dimensions. For a `GeoSeries`, the output is another series where
the child objects are sliced and the series dimensions index is now of the child 
dimensions combined. `slice` on a `GeoSeries` with no dimensions will slice along
the dimensions shared by both the series and child object.
"""
slice(x::Union{AbstractGeoArray,AbstractGeoStack}, dims) = slice(x, (dims,))
function slice(x::Union{AbstractGeoArray,AbstractGeoStack}, dims::Tuple)
    all(hasdim(x, dims)) || _errordimsnotfound(otherdims(dims, DD.dims(x)))
    seriesdims = DD.dims(x, dims)
    seriesdata = map(DD.dimwise_generators(seriesdims)) do ds
        view(x, ds...)
    end
    GeoSeries(seriesdata, seriesdims)
end
slice(ser::AbstractGeoSeries, dims) = cat(map(x -> slice(x, dims), ser)...; dims=dims)

@noinline _errordimsnotfound(dims) = 
    throw(ArgumentError("$(map(DD.dim2key, dims)) were dims not found in $(nameof(typeof(x)))"))
@noinline _errordimsnotfound(dims::Tuple{<:Any}) = 
    throw(ArgumentError("$(DD.dim2key(dims[1])) was dim not found in $(nameof(typeof(x)))"))

# Get the bounds wrapped in Dim(Between)
dimbounds(A::AbstractDimArray) = dimbounds(bounds, A)
function dimbounds(f::Function, A::AbstractDimArray)
    map(dims(A), f(A)) do dim, bounds
        # TODO is Between the right selector? do we need inclusive?
        basetypeof(dim)(Between(bounds))
    end
end

"""
    chunk(A::AbstractGeoArray)

Creat a GeoSeries of arrays matching the chunks of a chunked array. 

This may be useful for parallel or larger than memory applications.
"""
function chunk(A::AbstractGeoArray)
    gc = DiskArrays.eachchunk(A)
    ci = CartesianIndices(gc.chunkgridsize)
    data = collect(view(A, _chunk_inds(gc, I)...) for I in ci)
    GeoSeries(data, DD.basedims(dims(A)))
end

# See iterate(::GridChunks) in Diskarrays.jl
function _chunk_inds(g, ichunk) 
    outinds = map(ichunk.I, g.chunksize, g.parentsize, g.offset) do ic, cs, ps, of
        max((ic - 1) * cs + 1 -of, 1):min(ic * cs - of, ps)
    end
end
