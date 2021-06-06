using GeoData, Test, Statistics, Dates, CFTime, Plots
import ArchGDAL, NCDatasets
using GeoData: name, window, mode, span, sampling, val, Ordered, metadata, bounds,
               FileArray, FileStack, NCDfile
include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = maybedownload(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = maybedownload(joinpath(ncexamples, "test_echam_spectral.nc"))

stackkeys = (
    :abso4, :aclcac, :aclcov, :ahfcon, :ahfice, :ahfl, :ahfliac, :ahfllac,
    :ahflwac, :ahfres, :ahfs, :ahfsiac, :ahfslac, :ahfswac, :albedo, :albedo_nir,
    :albedo_nir_dif, :albedo_nir_dir, :albedo_vis, :albedo_vis_dif, :albedo_vis_dir,
    :alsobs, :alsoi, :alsol, :alsom, :alsow, :ameltdepth, :ameltfrac, :amlcorac,
    :ao3, :apmeb, :apmegl, :aprc, :aprl, :aprs, :aps, :az0i, :az0l, :az0w,
    :barefrac, :dew2, :drain, :evap, :evapiac, :evaplac, :evapwac, :fage, :friac,
    :geosp, :glac, :gld, :hyai, :hyam, :hybi, :hybm, :lsp, :q, :qres, :qvi, :relhum,
    :runoff, :sd, :seaice, :siced, :sicepdi, :sicepdw, :sicepres, :slm, :sn, :snacl,
    :snc, :sni, :snifrac, :snmel, :sofliac, :sofllac, :soflwac, :srad0, :srad0d,
    :srad0u, :sradl, :srads, :sradsu, :sraf0, :srafl, :srafs, :st, :svo, :t2max,
    :t2min, :temp2, :thvsig, :topmax, :tpot, :trad0, :tradl, :trads, :tradsu,
    :traf0, :trafl, :trafs, :trfliac, :trfllac, :trflwac, :tropo, :tsi, :tsicepdi,
    :tslm1, :tsurf, :tsw, :u10, :ustr, :ustri, :ustrl, :ustrw, :v10, :vdis, :vdisgw,
    :vstr, :vstri, :vstrl, :vstrw, :wimax, :wind10, :wl, :ws, :wsmx, :xi, :xivi,
    :xl, :xlvi
)

@testset "geoarray" begin
    ncarray = geoarray(ncsingle)

    @testset "open" begin
        @test all(open(A -> A[Y=1], ncarray) .=== ncarray[:, 1, :])
    end

    @testset "read" begin
        A = read(ncarray)
        @test A isa GeoArray
        @test parent(A) isa Array
    end

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test ncarray isa GeoArray
        @test val(dims(ncarray, Ti())) == DateTime360Day(2001, 1, 16):Month(1):DateTime360Day(2002, 12, 16)
        @test val(dims(ncarray, Y())) == -79.5:89.5
        @test val(dims(ncarray, X())) == 1.0:2:359.0
        @test bounds(ncarray) == (
            (0.0, 360.0), 
            (-80.0, 90.0), 
            (DateTime360Day(2001, 1, 1), DateTime360Day(2003, 1, 1)),
        )
    end

    @testset "dimensions" begin
        @test ndims(ncarray) == 3
        @test length.(dims(ncarray)) == (180, 170, 24)
        @test dims(ncarray) isa Tuple{<:X,<:Y,<:Ti}
        @test refdims(ncarray) == ()
        # TODO detect the time span, and make it Regular
        modes = (
             Mapped(Ordered(), Explicit(vcat((0.0:2.0:358.0)', (2.0:2.0:360.0)')), Intervals(Center()), EPSG(4326), EPSG(4326)),
             Mapped(Ordered(), Explicit(vcat((-80.0:89.0)', (-79.0:90.0)')), Intervals(Center()), EPSG(4326), EPSG(4326)),
             Sampled(Ordered(), Explicit(
                 vcat(permutedims(DateTime360Day(2001, 1, 1):Month(1):DateTime360Day(2002, 12, 1)), 
                      permutedims(DateTime360Day(2001, 2, 1):Month(1):DateTime360Day(2003, 1, 1)))
                ), Intervals(Center())
            )
        )
        @test val.(span(ncarray)) == val.(span.(modes))
        @test typeof(mode(ncarray)) == typeof(modes)
        @test bounds(ncarray) == ((0.0, 360.0), (-80.0, 90.0), (DateTime360Day(2001, 1, 1), DateTime360Day(2003, 1, 1)))
    end

    @testset "other fields" begin
        @test ismissing(missingval(ncarray))
        @test metadata(ncarray) isa Metadata{NCDfile}
        @test name(ncarray) == :tos
    end

    @testset "indexing" begin
        @test ncarray[Ti(1)] isa GeoArray{<:Any,2}
        @test ncarray[Y(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[X(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[X(1), Y(1), Ti(1)] isa Missing
        @test ncarray[X(30), Y(30), Ti(1)] isa Float32
        # Russia
        @test ncarray[X(50), Y(100), Ti(1)] isa Missing
        # Alaska
        @test ncarray[Y(Near(64.2008)), X(Near(149.4937)), Ti(1)] isa Missing
        @test ncarray[Ti(2), X(At(59.0)), Y(At(-50.5))] == ncarray[30, 30, 2] === 278.47168f0
    end

    @testset "indexing with reverse lat" begin
        if !haskey(ENV, "CI") # CI downloads fail. But run locally
            ncrevlat = maybedownload("ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst.v5/sst.mon.ltm.1981-2010.nc")
            ncrevlatarray = geoarray(ncrevlat; key=:sst, missingval=-9.96921f36)
            @test order(dims(ncrevlatarray, Y)) == Ordered(ReverseIndex(), ReverseArray(), ForwardRelation())
            @test ncrevlatarray[Y(At(40)), X(At(100)), Ti(1)] == missingval(ncrevlatarray)
            @test ncrevlatarray[Y(At(-40)), X(At(100)), Ti(1)] == ncrevlatarray[51, 65, 1] == 14.5916605f0
            @test val(span(ncrevlatarray, Ti)) == Month(1)
            @test val(span(ncrevlatarray, Ti)) isa Month # Not CompoundPeriod
        end
    end

    @testset "selectors" begin
        a = ncarray[Lon(At(21.0)), Lat(Between(50, 52)), Ti(Near(DateTime360Day(2002, 12)))]
        @test bounds(a) == ((50.0, 52.0),)
        x = ncarray[Lon(Near(150)), Lat(Near(30)), Ti(1)]
        @test x isa Float32
        dimz = Lon(Between(0.0, 360)), Lat(Between(-80, 90)), 
               Ti(Between(DateTime360Day(2001, 1, 1), DateTime360Day(2003, 01, 02)))
        @test size(ncarray[dimz...]) == (180, 170, 24)
        @test index(ncarray[dimz...]) == index(ncarray)
        nca = ncarray[Lat(Between(-80, -25)), Lon(Between(0, 180)), Ti(Contains(DateTime360Day(2002, 02, 20)))]
        @test size(nca) == (90, 55)
        @test index(nca, Lat) == index(ncarray[1:90, 1:55, 2], Lat)
        @test all(nca .=== ncarray[1:90, 1:55, 14])
    end

    @testset "conversion to GeoArray" begin
        geoA = ncarray[X(1:50), Y(20:20), Ti(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: Union{Missing,Float32}
        @time geoA isa GeoArray{Float32,1}
        @test dims(geoA) isa Tuple{<:X,<:Y}
        @test refdims(geoA) isa Tuple{<:Ti}
        @test metadata(geoA) == metadata(ncarray)
        @test ismissing(missingval(geoA))
        @test name(geoA) == :tos
    end

    @testset "save" begin
        @testset "to netcdf" begin
            # TODO save and load subset
            geoA = read(ncarray)
            metadata(geoA)
            @test size(geoA) == size(ncarray)
            filename = tempname() * ".nc"
            write(filename, geoA)
            saved = read(geoarray(filename))
            @test size(saved) == size(geoA)
            @test refdims(saved) == refdims(geoA)
            @test missingval(saved) === missingval(geoA)
            @test map(metadata.(dims(saved)), metadata.(dims(geoarray))) do s, g
                all(s .== g)
            end |> all
            @test_broken metadata(saved) == metadata(geoA)
            @test_broken all(metadata.(dims(saved)) .== metadata.(dims(geoA)))
            @test GeoData.name(saved) == GeoData.name(geoA)
            @test all(mode.(dims(saved)) .!= mode.(dims(geoA)))
            @test all(order.(dims(saved)) .== order.(dims(geoA)))
            @test all(typeof.(span.(dims(saved))) .== typeof.(span.(dims(geoA))))
            @test all(val.(span.(dims(saved))) .== val.(span.(dims(geoA))))
            @test all(sampling.(dims(saved)) .== sampling.(dims(geoA)))
            @test typeof(dims(saved)) <: typeof(dims(geoA))
            @test val(dims(saved)[3]) == val(dims(geoA)[3])
            @test all(val.(dims(saved)) .== val.(dims(geoA)))
            @test all(data(saved) .=== data(geoA))
            @test saved isa typeof(geoA)
            # TODO test crs
        end
        @testset "to gdal" begin
            gdalfilename = tempname() * ".tif"
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write(gdalfilename, nccleaned)
            gdalarray = geoarray(gdalfilename)
            # gdalarray WKT is missing one AUTHORITY
            # @test_broken crs(gdalarray) == convert(WellKnownText, EPSG(4326))
            # But the Proj representation is the same
            @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test bounds(gdalarray) == (bounds(nccleaned)..., (1, 1))
            # Tiff locus = Start, Netcdf locus = Center
            @test reverse(val(dims(gdalarray, Y))) .+ 0.5 ≈ val(dims(nccleaned, Y))
            @test val(dims(gdalarray, X)) .+ 1.0  ≈ val(dims(nccleaned, X))
            @test reverse(GeoArray(gdalarray); dims=Y()) ≈ nccleaned
        end
        @testset "to grd" begin
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write("testgrd.gri", nccleaned)
            grdarray = geoarray("testgrd.gri");
            @test crs(grdarray) == convert(ProjString, EPSG(4326))
            @test bounds(grdarray) == (bounds(nccleaned)..., (1, 1))
            @test val(dims(grdarray, Y)) ≈ val(dims(nccleaned, Y)) .- 0.5
            @test val(dims(grdarray, X)) ≈ val(dims(nccleaned, X)) .- 1.0
            @test GeoArray(grdarray) ≈ reverse(nccleaned; dims=Y)
        end
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), ncarray)
        # Test but don't lock this down too much
        @test occursin("GeoArray", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Time", sh)
    end

    @testset "plot" begin
        ncarray[Ti(1:3:12)] |> plot
        ncarray[Ti(1)] |> plot
        ncarray[Y(100), Ti(1)] |> plot
    end

end

@testset "Single file stack" begin
    ncstack = stack(ncmulti)

    @testset "load ncstack" begin
        @test ncstack isa GeoStack
        @test ismissing(missingval(ncstack))
        @test metadata(ncstack) isa Metadata{NCDfile}
        @test dims(ncstack[:abso4]) == dims(ncstack, (X, Y, Ti)) 
        @test refdims(ncstack) == ()
        # Loads child as a regular GeoArray
        @test_throws ErrorException ncstack[:not_a_key]
        @test ncstack[:albedo] isa GeoArray{<:Any,3}
        @test ncstack[:albedo, 2, 3, 1] isa Float32
        @test ncstack[:albedo, :, 3, 1] isa GeoArray{<:Any,1}
        @test dims(ncstack[:albedo]) isa Tuple{<:X,<:Y,<:Ti}
        @test keys(ncstack) isa NTuple{131,Symbol}
        @test keys(ncstack) == stackkeys
        @test first(keys(ncstack)) == :abso4
        @test metadata(ncstack) isa Metadata{NCDfile}
        @test metadata(ncstack)["institution"] == "Max-Planck-Institute for Meteorology"
        @test metadata(ncstack, :albedo) isa Metadata{NCDfile}
        @test metadata(ncstack, :albedo)["long_name"] == "surface albedo"
        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(ncstack[:albedo, Y(1:20)] , dims=Ti)) ==
              (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        geoA = ncstack[:albedo][Ti(4:6), X(1), Y(2)]
        @test geoA == ncstack[:albedo, Ti(4:6), X(1), Y(2)]
        @test size(geoA) == (3,)
    end

    @testset "read" begin
        st = read(ncstack)
        @test st isa GeoStack
        @test st.data isa NamedTuple
        @test first(st.data) isa Array
        GeoStack(st...)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = read(ncstack[:albedo]) .* 2
            copy!(geoA, ncstack, :albedo);
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoA == read(ncstack[:albedo])
        end
    end

    @testset "indexing" begin
        ncmultistack = stack(ncsingle)
        @test dims(ncmultistack[:tos]) isa Tuple{<:X,<:Y,<:Ti}
        @test ncmultistack[:tos] isa GeoArray{<:Any,3}
        @test ncmultistack[:tos, Ti(1)] isa GeoArray{<:Any,2}
        @test ncmultistack[:tos, Y(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncmultistack[:tos, 8, 30, 10] isa Float32
    end

    @testset "window" begin
        windowedstack = stack(ncmulti; window=(Y(1:5), X(1:5), Ti(1)))
        @test window(windowedstack) == (Y(1:5), X(1:5), Ti(1))
        windowedarray = windowedstack[:albedo]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] == 0.84936917f0
        windowedstack = stack(ncmulti; window=(Y(1:5), X(1:5), Ti(1:1)))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2, 1] == 0.84936917f0
        windowedstack = stack(ncmulti; window=(Ti(1),))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] ==  0.84936917f0
    end

    @testset "conversion to GeoStack" begin
        geostack = GeoStack(ncstack)
        @test Symbol.(Tuple(keys(geostack))) == keys(ncstack)
        smallstack = GeoStack(ncstack; keys=(:albedo, :evap, :runoff))
        @test keys(smallstack) == (:albedo, :evap, :runoff)
    end

    @testset "save" begin
        geostack = stack(ncstack);
        length(dims(geostack[:aclcac]))
        filename = tempname() * ".nc"
        write(filename, geostack);
        saved = GeoStack(stack(filename))
        @test keys(saved) == keys(geostack)
        @test metadata(saved)["advection"] == "Lin & Rood"
        @test metadata(saved) == metadata(geostack) == metadata(ncstack)
        @test all(first(DimensionalData.layers(saved)) .== first(DimensionalData.layers(geostack)))
    end

    @testset "show" begin
        ncstack = stack(ncmulti; window=(X(7:99), Y(3:90)));
        sh = sprint(show, MIME("text/plain"), ncstack)
        # Test but don't lock this down too much
        @test occursin("GeoStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Ti", sh)
        @test occursin(":tropo", sh)
        @test occursin(":tsurf", sh)
        @test occursin(":aclcac", sh)
        @test occursin("test_echam_spectral.nc", sh)
        @test occursin("window", sh)
        @test occursin("7:99", sh)
        @test occursin("3:90", sh)
    end

end

@testset "Multi file stack" begin
    ncstack = stack((tropo=ncmulti, tsurf=ncmulti, aclcac=ncmulti, albedo=ncmulti))
    @test length(ncstack) == 4
    @test dims(ncstack) isa Tuple{<:X,<:Y,<:Ti,<:Z}

    @testset "read" begin
        st = read(ncstack)
        @test st isa GeoStack
        @test st.data isa NamedTuple
        @test st.data[1] isa Array
        @test st.data[2] isa Array
    end

    @testset "child array properties" begin
        @test size(ncstack[:tropo]) == (192, 96, 8)
        @test ncstack[:tropo] isa GeoArray{Float32,3}
    end

    @testset "indexing" begin
        @test ncstack[:aclcac, Ti(1)] == ncstack[:aclcac][Ti(1)]
        @test typeof(ncstack[:aclcac, Ti(1)]) == typeof(ncstack[:aclcac][Ti(1)])
    end

    @testset "window" begin
        windowedstack = stack(ncmulti; window=(Y(1:5), X(1:5), Ti(1)))
        @test window(windowedstack) == (Y(1:5), X(1:5), Ti(1))
        windowedarray = windowedstack[:albedo]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] == 0.84936917f0
        windowedstack = stack(ncmulti; window=(Y(1:5), X(1:5), Ti(1:1)))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2, 1] == 0.84936917f0
        windowedstack = stack(ncmulti; window=(Ti(1),))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] ==  0.84936917f0
    end

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        geostack = GeoStack(ncstack)
        @test Symbol.(Tuple(keys(ncstack))) == keys(geostack)
        smallstack = GeoStack(ncstack; keys=(:tsurf,))
        @test keys(smallstack) == (:tsurf,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(GeoArray(ncstack[:tropo]))
            copy!(geoA, ncstack, :tropo)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test all(geoA .== GeoArray(ncstack[:tropo]))
        end
    end

    @testset "save" begin
        geoA = GeoArray(ncstack[:tsurf])
        filename = tempname() * ".nc"
        write(filename, ncstack)
        saved = read(geoarray(filename))
        @test_broken all(saved .== geoA)
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), ncstack)
        # Test but don't lock this down too much
        @test occursin("GeoStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Ti", sh)
        @test occursin(":tropo", sh)
        @test occursin(":tsurf", sh)
        @test occursin(":aclcac", sh)
    end

end

@testset "series" begin
    ncseries = series([ncmulti, ncmulti], (Ti,); child=stack)
    @testset "read" begin
        geoseries = read(ncseries)
        @test geoseries isa GeoSeries{<:GeoStack}
        @test geoseries.data isa Vector{<:GeoStack}
    end
    geoA = read(geoarray(ncmulti; key=:albedo))
    @test all(read(ncseries[Ti(1)][:albedo]) .== read(geoA))
    @test read(ncseries[Ti(1)][:albedo]) == read(geoA)
    @test all(read(ncseries[Ti(1)][:albedo]) .== read(geoA))
    @testset "modify" begin
        modified_series = modify(Array, ncseries)
        @test keys(modified_series) == keys(ncseries)
        @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{stackkeys,<:Tuple{<:Array{Float32,3},Vararg}}}}
    end
end

nothing
