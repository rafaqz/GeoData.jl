using GeoData, Test, Statistics, Dates, Plots
import NCDatasets, ArchGDAL
using GeoData: name, mode, window, bounds, FileArray, GRDfile, GDALfile

testpath = joinpath(dirname(pathof(GeoData)), "../test/")
include(joinpath(testpath, "test_utils.jl"))
const DD = DimensionalData

maybedownload("https://raw.githubusercontent.com/rspatial/raster/master/inst/external/rlogo.grd", "rlogo.grd")
maybedownload("https://github.com/rspatial/raster/raw/master/inst/external/rlogo.gri", "rlogo.gri")
stem = joinpath(testpath, "data/rlogo")
@test isfile(stem * ".grd")
@test isfile(stem * ".gri")
path = stem * ".gri"

@testset "Grd array" begin
    grdarray = geoarray(path)

    @testset "open" begin
        @test all(open(A -> A[Y=1], grdarray) .=== grdarray[:, 1, :])
        tempfile = tempname()
        cp(stem * ".grd", tempfile * ".grd")
        cp(stem * ".gri", tempfile * ".gri")
        grdwritearray = geoarray(tempfile * ".gri")
        open(grdwritearray; write=true) do A
            A .*= 2
        end
        @test geoarray(tempfile * ".gri") == grdarray .* 2
    end

    @testset "read" begin
        A = read(grdarray)
        @test A isa GeoArray
        @test parent(A) isa Array
    end

    @testset "array properties" begin
        @test grdarray isa GeoArray{Float32,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(grdarray), X))) == 101
        @test ndims(grdarray) == 3
        @test dims(grdarray) isa Tuple{<:X,<:Y,<:Band}
        @test refdims(grdarray) == ()
        @test bounds(grdarray) == ((0.0, 101.0), (0.0, 77.0), (1, 3))
    end

    @testset "other fields" begin
        @test missingval(grdarray) == -3.4f38
        @test metadata(grdarray) isa Metadata{GRDfile}
        @test name(grdarray) == Symbol("red:green:blue")
        @test label(grdarray) == "red:green:blue"
        @test units(grdarray) == nothing
        customgrdarray = geoarray(path; name=:test, mappedcrs=EPSG(4326));
        @test name(customgrdarray) == :test
        @test label(customgrdarray) == "test"
        @test mappedcrs(dims(customgrdarray, Y)) == EPSG(4326)
        @test mappedcrs(dims(customgrdarray, X)) == EPSG(4326)
        @test mappedcrs(customgrdarray) == EPSG(4326)
        proj = ProjString("+proj=merc +datum=WGS84")
        @test crs(dims(customgrdarray, Y)) == proj
        @test crs(dims(customgrdarray, X)) == proj
        @test crs(customgrdarray) == proj
    end

    @testset "getindex" begin
        @test grdarray[Band(1)] isa GeoArray{Float32,2}
        @test grdarray[Y(1), Band(1)] isa GeoArray{Float32,1}
        @test grdarray[X(1), Band(1)] isa GeoArray{Float32,1}
        @test grdarray[X(50), Y(30), Band(1)] == 115.0f0
        @test grdarray[1, 1, 1] == 255.0f0
        @test grdarray[Y(At(20)), X(At(20)), Band(3)] == 255.0f0
        @test grdarray[Y(Contains(60)), X(Contains(20)), Band(1)] == 255.0f0
    end

    @testset "selectors" begin
        geoA = grdarray[Y(Contains(3)), X(:), Band(1)]
        @test geoA isa GeoArray{Float32,1}
        @test grdarray[X(Contains(20)), Y(Contains(10)), Band(1)] isa Float32
    end

    @testset "conversion to GeoArray" begin
        geoA = grdarray[X(1:50), Y(1:1), Band(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: Float32
        @time geoA isa GeoArray{Float32,1}
        @test dims(geoA) isa Tuple{<:X,Y}
        @test refdims(geoA) isa Tuple{<:Band}
        @test metadata(geoA) == metadata(grdarray)
        @test missingval(geoA) == -3.4f38
        @test name(geoA) == Symbol("red:green:blue")
    end

    @testset "save" begin
        @testset "2d" begin
            filename2 = tempname() * ".gri"
            write(filename2, grdarray[Band(1)])
            saved = read(geoarray(filename2))
            # 1 band is added again on save
            @test size(saved) == size(grdarray[Band(1:1)])
            @test data(saved) == data(grdarray[Band(1:1)])
        end

        @testset "3d with subset" begin
            geoA = grdarray[1:100, 1:50, 1:2]
            filename = tempname() * ".grd"
            write(filename, GRDfile, geoA)
            saved = read(geoarray(filename))
            @test size(saved) == size(geoA)
            @test refdims(saved) == ()
            @test bounds(saved) == bounds(geoA)
            @test size(saved) == size(geoA)
            @test missingval(saved) === missingval(geoA)
            @test metadata(saved) != metadata(geoA)
            @test metadata(saved)["creator"] == "GeoData.jl"
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoA)))
            @test name(saved) == name(geoA)
            @test all(mode.(dims(saved)) .== mode.(dims(geoA)))
            @test dims(saved) isa typeof(dims(geoA))
            @test all(val.(dims(saved)) .== val.(dims(geoA)))
            @test all(mode.(dims(saved)) .== mode.(dims(geoA)))
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoA)))
            @test dims(saved) == dims(geoA)
            @test all(data(saved) .=== data(geoA))
            @test saved isa typeof(geoA)
            @test data(saved) == data(geoA)
        end

        @testset "to netcdf" begin
            filename2 = tempname() * ".nc"
            span(grdarray[Band(1)])
            write(filename2, grdarray[Band(1)])
            saved = read(geoarray(filename2; crs=crs(grdarray)))
            @test size(saved) == size(grdarray[Band(1)])
            @test replace_missing(saved, missingval(grdarray)) ≈ reverse(grdarray[Band(1)]; dims=Y)
            @test replace_missing(saved, missingval(grdarray)) ≈ reverse(grdarray[Band(1)]; dims=Y)
            @test index(saved, X) ≈ index(grdarray, X) .+ 0.5
            @test index(saved, Y) ≈ index(grdarray, Y) .+ 0.5
            @test bounds(saved, Y) == bounds(grdarray, Y)
            @test bounds(saved, X) == bounds(grdarray, X)
        end

        @testset "to gdal" begin
            # No Band
            gdalfilename = tempname() * ".tif"
            write(gdalfilename, GDALfile, grdarray[Band(1)])
            gdalarray = geoarray(gdalfilename)
            # @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test val(dims(gdalarray, X)) ≈ val(dims(grdarray, X))
            @test reverse(val(dims(gdalarray, Y))) ≈ val(dims(grdarray, Y))
            @test GeoArray(gdalarray) ≈ permutedims(grdarray[Band(1)], [X(), Y()])
            # 3 Bands
            gdalfilename2 = tempname() * ".tif"
            write(gdalfilename2, grdarray)
            gdalarray2 = geoarray(gdalfilename2)
            @test all(GeoArray(gdalarray2) .== GeoArray(grdarray))
            @test val(dims(gdalarray2, Band)) == 1:3
        end

    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), grdarray)
        # Test but don't lock this down too much
        @test occursin("GeoArray", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
    end

    @testset "plot" begin
        grdarray |> plot
        grdarray[Band(1)] |> plot
    end

end

@testset "Grd stack" begin
    grdstack = stack((a=path, b=path))

    @test length(grdstack) == 2
    @test dims(grdstack) isa Tuple{<:X,<:Y,<:Band}

    @testset "read" begin
        st = read(grdstack)
        @test st isa GeoStack
        @test st.data isa NamedTuple
        @test first(st.data) isa Array
    end

    @testset "indexing" begin
        @test grdstack[:a][Y(20), X(20), Band(3)] == 70.0f0
        @test grdstack[:a][Y([2,3]), X(40), Band(2)] == [240.0f0, 246.0f0]
    end

    @testset "child array properties" begin
        @test size(grdstack[:a]) == size(GeoArray(grdstack[:a])) == (101, 77, 3)
        @test grdstack[:a] isa GeoArray{Float32,3}
    end

    @testset "window" begin
        windowedstack = stack((a=path, b=path); window=(Y(1:5), X(1:5), Band(1)))
        @test window(windowedstack) == (Y(1:5), X(1:5), Band(1))
        windowedarray = windowedstack[:a]
        @test windowedarray isa GeoArray{Float32,2}
        @test length.(dims(windowedarray)) == (5, 5)
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[1, 2] == 255.0f0
        windowedstack = stack((a=path, b=path); window=(Y(1:5), X(1:5), Band(1:1)))
        windowedarray = windowedstack[:b]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[1, 2, 1] == 255.0f0
        windowedstack = stack((a=path, b=path); window=(Band(1),));
        windowedarray = windowedstack[:b]
        @test windowedarray[1:3, 2:2] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[30, 30] == 185.0f0
        windowedarray |> plot
    end

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        geostack = GeoStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(geostack)
        smallstack = GeoStack(grdstack; keys=(:a,))
        @test keys(smallstack) == (:a,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(GeoArray(grdstack[:a]))
            copy!(geoA, grdstack, :a)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoA == GeoArray(grdstack[:a])
        end
    end

    @testset "save" begin
        geoA = GeoArray(grdstack[:b])
        filename = tempname() * ".grd"
        write(filename, grdstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = read(geoarray(filename_b))
        @test typeof(saved) == typeof(geoA)
        @test data(saved) == data(geoA)
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), grdstack)
        # Test but don't lock this down too much
        @test occursin("GeoStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
        @test occursin(":a", sh)
        @test occursin(":b", sh)
    end

end

@testset "Grd series" begin
    grdseries = series([path, path], (Ti,); mappedcrs=EPSG(4326))
    @test grdseries[Ti(1)] == geoarray(path; mappedcrs=EPSG(4326))
    stacks = [stack((a=path, b=path); mappedcrs=EPSG(4326))]

    grdseries2 = GeoSeries(stacks, (Ti,))
    @test all(grdseries2[Ti(1)][:a] .== geoarray(path; mappedcrs=EPSG(4326), name=:test))
    modified_ser = modify(x -> Array(1.0f0x), grdseries2)
    @test typeof(modified_ser) <: GeoSeries{<:GeoStack{<:NamedTuple{(:a,:b),<:Tuple{<:Array{Float32,3},Vararg}}},1}

    @testset "read" begin
        geoseries = read(grdseries2)
        @test geoseries isa GeoSeries{<:GeoStack}
        @test geoseries.data isa Vector{<:GeoStack}
        @test geoseries.data isa Vector{<:GeoStack}
        @test first(geoseries.data[1].data) isa Array 
    end
end


nothing
