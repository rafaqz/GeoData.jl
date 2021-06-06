export reproject, convertmode, mappedindex, mappedbounds

# These methods are only available if ArchGDAL is loaded.
# Otherwise Projected selector crs field must be `nothing`,
# and no reprojection can occur.

function DD._sel2indices(mode::Projected, dim::Dimension, sel::Contains)
    selval = reproject(mappedcrs(mode), crs(mode), dim, val(sel))
    DD.contains(dim, rebuild(sel, selval))
end
function DD._sel2indices(mode::Projected, dim::Dimension, sel::At)
    selval = reproject(mappedcrs(mode), crs(mode), dim, val(sel))
    DD.at(dim, rebuild(sel, selval))
end
function DD._sel2indices(mode::Projected, dim::Dimension, sel::Between)
    selval = map(v -> reproject(mappedcrs(mode), crs(mode), dim, v), val(sel))
    DD.between(dim, rebuild(sel, selval))
end

"""
    reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a single dimension at a time.
"""
reproject(source, target, dim::Dimension, val) = val
function reproject(source::GeoFormat, target::GeoFormat, dim::X, val::Number)
    AG.reproject((val, zero(val)), source, target; order=:trad)[1]
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, val::Number)
    AG.reproject((zero(val), val), source, target; order=:trad)[2]
end
function reproject(source::GeoFormat, target::GeoFormat, ::X, vals::AbstractArray)
    rep = AG.reproject(map(v -> (v, zero(v)), vals), source, target; order=:trad)
    map(r -> r[1], rep)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, vals::AbstractArray)
    rep = AG.reproject(map(v -> (zero(v), v), vals), source, target; order=:trad)
    map(r -> r[2], rep)
end
function reproject(source::GeoFormat, target::GeoFormat, ::X, vals::Tuple)
    reps = AG.reproject([(v, zero(v)) for v in vals], source, target; order=:trad)
    Tuple(r[1] for r in reps)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, vals::Tuple)
    reps = AG.reproject([(zero(v), v) for v in vals], source, target; order=:trad)
    Tuple(r[2] for r in reps)
end

function convertmode(dstmode::Type{Mapped}, srcmode::Type{Projected}, dim::Dimension)
    m = mode(dim)
    newindex = reproject(crs(m), mappedcrs(m), dim, index(dim))
    newbounds = reproject(crs(m), mappedcrs(m), dim, DD.dim2boundsmatrix(dim))
    newmode = Mapped(
        order=order(m),
        span=Explicit(newbounds),
        sampling=sampling(m),
        crs=crs(m),
        mappedcrs=mappedcrs(m),
    )
    rebuild(dim; val=newindex, mode=newmode)
end
function convertmode(dstmode::Type{Projected}, srcmode::Type{Mapped}, dim::Dimension)
    m = mode(dim)
    newindex = _projectedrange(m, dim)
    newmode = Projected(
        order=order(m),
        span=Regular(step(newindex)), sampling=sampling(m),
        crs=crs(m),
        mappedcrs=mappedcrs(m),
    )
    rebuild(dim; val=newindex, mode=newmode)
end

_projectedrange(::Projected, dim) = LinRange(first(dim), last(dim), length(dim))
_projectedrange(m::Mapped, dim) = _projectedrange(span(m), crs(m), m, dim)
_projectedrange(span, crs, m::Mapped, dim) = begin
    start, stop = reproject(mappedcrs(m), crs, dim, [first(dim), last(dim)])
    LinRange(start, stop, length(dim))
end
_projectedrange(::Regular, crs::Nothing, ::Mapped, dim) =
    LinRange(first(dim), last(dim), length(dim))
_projectedrange(::T, crs::Nothing, ::Mapped, dim) where T<:Union{Irregular,Explicit} =
    error("Cannot convert a Mapped $T index to Projected when crs is nothing")


projectedbounds(mode::Mapped, dim) = projectedbounds(crs(mode), mode, dim)
projectedbounds(crs::Nothing, mode::Mapped, dim) =
    error("No projection crs attached to $(name(dim)) dimension")
projectedbounds(crs::GeoFormat, mode::Mapped, dim) =
    reproject(mappedcrs(mode), crs, dim, bounds(dim))

projectedindex(mode::Mapped, dim) = projectedindex(crs(mode), mode, dim)
projectedindex(crs::Nothing, mode::Mapped, dim) =
    error("No projection crs attached to $(name(dim)) dimension")
projectedindex(crs::GeoFormat, mode::Mapped, dim) =
    reproject(mappedcrs(dim), crs, dim, index(dim))

mappedbounds(mode::Projected, dim) = mappedbounds(mappedcrs(mode), mode, dim)
mappedbounds(mappedcrs::Nothing, mode::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
mappedbounds(mappedcrs::GeoFormat, mode::Projected, dim) =
    reproject(crs(mode), mappedcrs, dim, bounds(dim))

mappedindex(mode::Projected, dim) = mappedindex(mappedcrs(mode), mode, dim)
mappedindex(mappedcrs::Nothing, mode::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
mappedindex(mappedcrs::GeoFormat, mode::Projected, dim) =
    reproject(crs(dim), mappedcrs, dim, index(dim))
