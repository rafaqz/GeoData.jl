using GeoData, Test, Aqua, SafeTestsets

if VERSION >= v"1.5.0"
    # Aqua.test_ambiguities([GeoData, Base, Core])
    # Aqua.test_unbound_args(GeoData)
    # Aqua.test_stale_deps(GeoData)
    Aqua.test_undefined_exports(GeoData)
    Aqua.test_project_extras(GeoData)
    Aqua.test_deps_compat(GeoData)
    Aqua.test_project_toml_formatting(GeoData)
end

@time @safetestset "array" begin include("array.jl") end
@time @safetestset "stack" begin include("stack.jl") end
@time @safetestset "series" begin include("series.jl") end
@time @safetestset "utils" begin include("utils.jl") end
@time @safetestset "set" begin include("set.jl") end
@time @safetestset "reproject" begin include("reproject.jl") end
@time @safetestset "aggregate" begin include("aggregate.jl") end
@time @safetestset "methods" begin include("methods.jl") end
@time @safetestset "resample" begin include("resample.jl") end
# Only test SMAP locally for now
if !haskey(ENV, "CI")
    @time @safetestset "smap" begin include("sources/smap.jl") end
end
if !Sys.iswindows()
    # GDAL Environment vars need to be set manually for windows, so skip for now
    @time @safetestset "gdal" begin include("sources/gdal.jl") end
    @time @safetestset "grd" begin include("sources/grd.jl") end
end
@time @safetestset "ncdatasets" begin include("sources/ncdatasets.jl") end
@time @safetestset "rasterdatasources" begin include("sources/rasterdatasources.jl") end
@time @safetestset "run examples" begin include("../examples/examples.jl") end
