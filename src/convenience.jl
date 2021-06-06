"""
    geoarray(filename; kw...) => AbstractGeoArray

Load a file path as `AbstractGeoArray`.

# Keywords

Passed to the constructor for the file type, and commmonly include:

- `crs`: crs to use instead of the detected crs.
- `mappedcrs`: CRS format like `EPSG(4326)` used in `Selectors` like `Between` and `At`, and
    for plotting. Can be any CRS `GeoFormat` from GeoFormatTypes.jl, like `WellKnownText`.
- `name`: `Symbol` name for the array.
- `dims`: `Tuple` of `Dimension`s for the array. Detected automatically, but can be passed in.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from.
- `missingval`: Value reprsenting missing values. Detected automatically when possible, but
    can be passed it.
- `metadata`: `Metadata` object for the array. Detected automatically but can be passed in.
"""
function geoarray end

geoarray(filename::AbstractString; kw...) = GeoArray(filename; kw...)

"""
    Base.write(filename::AbstractString, A::AbstractGeoArray; kw...)

Write any [`AbstractGeoArray`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.
"""
function Base.write(
    filename::AbstractString, A::AbstractGeoArray; kw...
)
    write(filename, _sourcetype(filename), A; kw...)
end

"""
    stack(filename; kw...) => AbstractGeoStack
    stack(filenames::NamedTuple; kw...) => AbstractGeoStack

Load a file path or a `NamedTuple` of paths as an `AbstractGeoStack`.

# Arguments

- `filename`: A `NamedTuple` of stack keys and `String` filenames.

# Keywords

Passed to the constructor for the file type, and commmonly include:

- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
    contained arrays when they are accessed.
- `metadata`: `Metadata` as object.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = stack(files; mappedcrs=EPSG(4326))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
function stack end

stack(args...; kw...) = GeoStack(args...; kw...)

"""
    Base.write(filename::AbstractString, T::Type{<:AbstractGeoArray}, s::AbstractGeoStack)

Write any [`AbstractGeoStack`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.

If the source can't be saved as a stack-like object, individual array layers will be saved.
"""
function Base.write(filename::AbstractString, s::AbstractGeoStack; kw...)
    base, ext = splitext(filename)
    T = _sourcetype(filename)
    if cansavestack(T)
        write(filename, _sourcetype(filename), s; kw...)
    else
        # Otherwise write separate arrays
        for key in keys(s)
            fn = joinpath(string(base, "_", key, ext))
            write(fn, _sourcetype(filename), s[key])
        end
    end
end

cansavestack(T) = false

"""
    series(dirpath::AbstractString, dims; ext, child=geoarray, kw...) => AbstractGeoSeries
    series(filenames::AbstractVector{String}, dims; child=geoarray, kw...) => AbstractGeoSeries

Load a whole folder vector of filepaths as a `AbstractGeoSeries`.
`kw` are passed to the constructor.

# Arguments

- `dims`: A tuple of dimensions, possibly holding an index matching the files in the directory.
    By default it will be a numbered `Dim{:series}`.

# Keywords

- `child`: function to load the child object - may be `geoarray` or `stack`,
    `geoarray` by default.
- `ext`: an extension for the files in a directory. If empty, all files are loaded,
    and should be the same type.
"""
function series end

function series(dirpath::AbstractString, dims=(Dim{:series}(),); ext=nothing, child=geoarray, kw...)
    filepaths = filter_ext(dirpath, ext)
    series(filepaths, dims; child=child, kw...)
end
function series(filepaths::AbstractVector{<:AbstractString}, dims=(Dim{:series}(),); child=geoarray, kw...)
    GeoSeries(filepaths, dims; child=child, kw...)
end

# Support methods

const EXT = Dict(GRDfile=>(".grd", ".gri"), NCDfile=>(".nc",), SMAPfile=>(".h5",))
const REV_EXT = Dict(".grd"=>GRDfile, ".gri"=>GRDfile, ".nc"=>NCDfile, ".h5"=>SMAPfile)

_sourcetype(filename::AbstractString) = get(REV_EXT, splitext(filename)[2], GDALfile)

function _read(f, filename::AbstractString; kw...)
    ext = splitext(filename)[2]
    source = get(REV_EXT, ext, GDALfile)
    _read(f, source, filename; kw...)
end
