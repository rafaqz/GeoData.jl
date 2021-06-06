struct GeoPlot end
struct GeoZPlot end

# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractGeoArray)
    ddplot(A) = DimArray(A; dims=_maybe_mapped(dims(A)))
    A = GeoArray(A)
    if !(get(plotattributes, :seriestype, :none) in (:none, :heatmap))
        DD.DimensionalPlot(), ddplot(A)
    elseif all(hasdim(A, (SpatialDim, SpatialDim)))
        # Heatmap or multiple heatmaps. Use GD recipes.
        GeoPlot(), _prepare(A)
    elseif hasdim(A, ZDim) && ndims(A) == 1
        # Z dim plot, but for spatial data we want Z on the Y axis
        GeoZPlot(), _prepare(A)
    else
        DD.DimensionalPlot(), ddplot(A)
    end
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::GeoPlot, A::GeoArray{T,3,<:Tuple{<:SpatialDim,<:SpatialDim,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :title --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                subplot := i
                GeoPlot(), A[:, :, i]
                slice = A[:, :, i]
                ys, xs = map(_prepare, dims(A))
                xs, ys, parent(slice)
            end
        end
    else
        GeoPlot(), A[:, :, 1]
    end
end
# Plot a sinlge 2d map
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:SpatialDim,<:SpatialDim}}) where T
    # If colorbar is close to symmetric (< 25% difference) use a symmetric 
    # colormap and set symmetric limits so zero shows up as a neutral color.
    A_min, A_max = extrema(skipmissing(A))
    if (A_min + A_max) / abs(A_max - A_min) < 0.25
        A_limit = max(abs(A_min), abs(A_max))
        clims = (-A_limit, A_limit)
        :seriescolor --> :balance
    else
        clims = A_min, A_max
    end

    yguide, xguide = label(dims(A))

    :seriestype --> :heatmap
    :title --> "$(_maybename(A)) $(DD._refdims_title(A))"
    :xguide --> xguide
    :yguide --> yguide
    :clims --> clims
    :colorbar_title --> name(A)

    ys, xs = map(_prepare, dims(A))

    if get(plotattributes, :seriestype, :none) == :contourf
        :linewidth --> 0
        :levels --> range(clims[1], clims[2], length=20)
        xs, ys, clamp.(A, clims[1], clims[2])
    else
        xs, ys, parent(A)
    end
end

# # Plot a vertical 1d line
@recipe function f(::GeoZPlot, A::GeoArray{T,1}) where T
    z_dim = dims(A, ZDim)
    yguide = label(z_dim)
    xguide = label(A)
    :title --> "$(_maybename(A))$(DD._refdims_title(A))"
    :xguide --> xguide
    :yguide --> yguide
    :label --> ""
    z = map(_prepare, dims(A))
    parent(A), z
end

# Plots heatmaps pixels are centered.
# So we should center, and use the projected value.
_prepare(d::Dimension) = d |> _maybe_shift |> _maybe_mapped |> index
# Convert arrays to a consistent missing value and Forward array order
_prepare(A::AbstractGeoArray) =
    replace_missing(A, missing) |>
    A -> reorder(A, ForwardIndex) |>
    A -> reorder(A, ForwardRelation) |>
    A -> permutedims(A, DD.commondims(>:, (ZDim, YDim, XDim, TimeDim, Dimension), dims(A)))

_maybename(A) = _maybename(name(A))
_maybename(n::Symbol) = n == Symbol("") ? "" : string(n, ": ")
_maybename(n::Name{N}) where N = _maybename(N) 
_maybename(n::NoName) = ""

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, eltype(A)(NaN))
_maybe_replace_missing(A) = A

_maybe_shift(d) = _maybe_shift(sampling(d), d)
_maybe_shift(::Intervals, d) = DD.maybeshiftlocus(Center(), d)
_maybe_shift(sampling, d) = d

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(mode(dim), dim)
_maybe_mapped(mode::IndexMode, dim::Dimension) = dim
_maybe_mapped(mode::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(mode), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertmode(Mapped, dim)
