
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{<:Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{<:Int}}}

"""
    aggregate(method, object, scale)

Aggregate array, or all arrays in a stack or series, by some scale.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `object`: Object to aggregate, like `AbstractGeoSeries`, `AbstractStack`,
  `AbstractGeoArray`, `Dimension`
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.
"""
function aggregate end
"""
    aggregate(method, series::AbstractGeoSeries, scale)

Aggregate an [`AbstractGeoSeries`](@ref) by `scale` using `method`.

Returns a [`GeoSeries`](ref).
"""
function aggregate(
    method, series::AbstractGeoSeries, scale, args...; progress=true, kwargs...
)
    f = i -> aggregate(method, series[i], scale, args...; progress=false, kwargs...)
    data = if progress
        @showprogress "Aggregating series..." map(f, 1:length(series))
    else
        map(f, 1:length(series))
    end
    return rebuild(series, data)
end
"""
    aggregate(method, stack::AbstractGeoStack, scale)

Aggregate an [`AbstractGeoStack`](@ref) by `scale` using `method`.

Returns a [`GeoStack`](ref).
"""
function aggregate(
    method, stack::AbstractGeoStack, scale; keys=keys(stack), progress=true
)
    f = key -> aggregate(method, stack[key], scale)
    keys_nt = NamedTuple{keys}(keys)
    arrays = if progress
        @showprogress "Aggregating stack..." map(f, keys_nt)
    else
        map(f, keys_nt)
    end
    return GeoStack(arrays)
end
# DimensionalData methods
"""
    aggregate(method, src::AbstractDimArray, scale)

Aggregate an `AbstractDimArray` by `scale` using `method`.
"""
function aggregate(method, src::AbstractDimArray, scale)
    aggregate!(method, alloc_ag(method, src, scale), src, scale)
end
"""
    aggregate(method, dim::Dimension, scale)

Aggregate a `Dimension` by `scale` using `method`.
"""
function aggregate(method, dim::Dimension, scale)
    intscale = _scale2int(dim, scale)
    intscale == 1 && return dim
    start, stop = _endpoints(dim, method, intscale)
    return rebuild(dim, val(dim)[start:scale:stop], aggregate(method, mode(dim), intscale))
end
"""
    aggregate(method, dim::IndexMode, scale)

Aggregate an `IndexMode` by `scale` using `method`.
"""
aggregate(method, mode::IndexMode, scale) = mode
function aggregate(method, mode::AbstractSampled, scale)
    rebuild(mode; span=aggregate(method, span(mode), scale))
end
"""
    aggregate(method, dim::Span, scale)

Aggregate a `Span` by `scale` using `method`.
"""
aggregate(method, span::Span, scale) = span
aggregate(method, span::Regular, scale) = Regular(val(span) * scale)
"""
    aggregate!(method, dst::AbstractDimArray, src::AbstractDimArray, scale)

Aggregate array `src` to array `dst` by `scale`, using `method`.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.
"""
function aggregate!(locus::Locus, dst::AbstractDimArray, src, scale)
    aggregate!((locus,), dst, src, scale)
end
function aggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractDimArray, src, scale)
    intscale = _scale2int(dims(src), scale)
    offsets = _agoffset.(loci, intscale)
    for I in CartesianIndices(dst)
        dst[I] = src[(upsample.(Tuple(I), intscale) .+ offsets)...]
    end
    return dst
end
# Function/functor methods
function aggregate!(f, dst::AbstractDimArray, src, scale)
    intscale = _scale2int(dims(src), scale)
    for I in CartesianIndices(dst)
        topleft = upsample.(Tuple(I), intscale)
        bottomright = topleft .+ intscale .- 1
        dst[I] = f(view(src, map(:, topleft, bottomright)...))
    end
    return dst
end


"""
    disaggregate(method, object, scale)

Disaggregate array, or all arrays in a stack or series, by some scale.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `object`: Object to aggregate, like `AbstractGeoSeries`, `AbstractStack`,
  `AbstractGeoArray`, `Dimension`
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.
"""
function disaggregate end
"""
    disaggregate(method, series::AbstractGeoSeries, scale)

Disaggregate an [`AbstractGeoSeries`](@ref) by `scale` using `method`.
"""
function disaggregate(method, series::AbstractGeoSeries, scale; progress=true, kwargs...)
    f = i -> disaggregate(method, series[i], scale; progress=false, kwargs...)
    return if progress
        @showprogress "Disaggregating series..." map(f, 1:length(series))
    else
        map(f, 1:length(series))
    end
end
"""
    disaggregate(method, stack::AbstractGeoStack, scale)

Disaggregate an [`AbstractGeoStack`](@ref) by `scale` using `method`.
"""
function disaggregate(method, stack::AbstractGeoStack, scale;
    keys=keys(stack), progress=true
)
    f = key -> disaggregate(method, stack[key], scale)
    keys_nt = NamedTuple{keys}(keys)
    arrays = if progress
        @showprogress "Disaggregating stack..." map(f, keys_nt)
    else
        map(f, keys_nt)
    end
    return GeoStack(arrays)
end
# DimensionalData methods
"""
    disaggregate(method, src::AbstractDimArray, scale)

Disaggregate an `AbstractDimArray` by `scale` using `method`.
"""
function disaggregate(method, src::AbstractDimArray, scale)
    disaggregate!(method, alloc_disag(method, src, scale), src, scale)
end
"""
    disaggregate(method, dim::Dimension, scale)

Disaggregate a `Dimension` by `scale` using `method`.
"""
function disaggregate(locus::Locus, dim::Dimension, scale)
    intscale = _scale2int(dim, scale)
    intscale == 1 && return dim
    len = length(dim) * intscale
    step_ = step(mode(dim)) / intscale
    start = dim[1] - _agoffset(locus, intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    return rebuild(dim, index, disaggregate(locus, mode(dim), intscale))
end
"""
    disaggregate(method, dim::IndexMode, scale)

Disaggregate an `IndexMode` by `scale` using `method`.
"""
disaggregate(method, mode::IndexMode, scale) = mode
disaggregate(method, mode::AbstractSampled, scale) =
    rebuild(mode; span=disaggregate(method, span(mode), scale))
"""
    disaggregate(method, dim::Span, scale)

Disaggregate a `Span` by `scale` using `method`.
"""
disaggregate(method, span::Span, scale) = span
disaggregate(method, span::Regular, scale) = Regular(val(span) / scale)
"""
    disaggregate!(method, dst::AbstractDimArray, src::AbstractDimArray, scale)

Disaggregate array `src` to array `dst` by some scale, using `method`.

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.
"""
function disaggregate!(locus::Locus, dst::AbstractDimArray, src, scale)
    disaggregate!((locus,), dst, src, scale)
end
function disaggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractDimArray, src, scale)
    intscale = _scale2int(dims(src), scale)
    for I in CartesianIndices(dst)
        dst[I] = src[(downsample.(Tuple(I), intscale))...]
    end
    return dst
end

# Allocate an array of the correct size to aggregate `A` by `scale`
alloc_ag(method, A::AbstractDimArray, scale) = alloc_ag((method,), A, scale)
function alloc_ag(method::Tuple, A::AbstractDimArray, scale)
    intscale = _scale2int(dims(A), scale)
    # Aggregate the dimensions
    dims_ = aggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    data_ = similar(data(A), map(length, dims_)...)
    return rebuild(A; data=data_, dims=dims_)
end

# Allocate an array of the correct size to disaggregate `A` by `scale`
alloc_disag(method, A::AbstractDimArray, scale) = alloc_disag((method,), A, scale)
function alloc_disag(method::Tuple, A::AbstractDimArray, scale)
    intscale = _scale2int(dims(A), scale)
    dims_ = disaggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    data_ = similar(data(A), map(length, dims_)...)
    return rebuild(A; data=data_, dims=dims_)
end


# Convert indicies from the aggregated array to the larger original array.
upsample(index::Int, scale::Int) = (index - 1) * scale + 1
upsample(index::Int, scale::Colon) = index

# Convert indicies from the original array to the aggregated array.
downsample(index::Int, scale::Int) = (index - 1) ÷ scale + 1
downsample(index::Int, scale::Colon) = index

# Convert scale or tuple of scale to integer using dims2indices
_scale2int(dims, scale::Tuple) = DD.dims2indices(dims, scale)
_scale2int(dims, scale::Int) = scale
_scale2int(dims, scale::Colon) = 1

_agoffset(locus::Locus, dim::Dimension, scale) = _agoffset(locus, scale)
_agoffset(method, dim::Dimension, scale) = _agoffset(locus(dim), scale)
_agoffset(locus::Start, scale) = 0
_agoffset(locus::End, scale) = scale - 1
_agoffset(locus::Center, scale) = scale ÷ 2

function _endpoints(dim, method, scale)
    start = firstindex(dim) + _agoffset(method, dim, scale)
    stop = (length(dim) ÷ scale) * scale
    return start, stop
end
