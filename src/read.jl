"""
    read(A::AbstractGeoArray)
    read(A::AbstractGeoStack)
    read(A::AbstractGeoSeries)

`read` will move a GeoData.jl object completely to memory, or make
a copy if already in memory.
```
"""
Base.read(x::Union{AbstractGeoArray,AbstractGeoStack,AbstractGeoSeries}) = modify(Array, x)
