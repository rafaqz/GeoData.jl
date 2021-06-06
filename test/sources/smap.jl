using GeoData, Test, Statistics, Dates, Plots
import ArchGDAL, NCDatasets, HDF5, CFTime
using GeoData: Time, window, name, layerkeys, SMAPfile

testpath = joinpath(dirname(pathof(GeoData)), "../test/")
include(joinpath(testpath, "test_utils.jl"))

# TODO example files without a login requirement
path1 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5")
path2 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5")

smapkeys = (
    :baseflow_flux, :heat_flux_ground, :heat_flux_latent, :heat_flux_sensible, 
    :height_lowatmmodlay, :land_evapotranspiration_flux, :land_fraction_saturated, 
    :land_fraction_snow_covered, :land_fraction_unsaturated, :land_fraction_wilting, 
    :leaf_area_index, :net_downward_longwave_flux, :net_downward_shortwave_flux, 
    :overland_runoff_flux, :precipitation_total_surface_flux, :radiation_longwave_absorbed_flux, 
    :radiation_shortwave_downward_flux, :sm_profile, :sm_profile_pctl, :sm_profile_wetness, 
    :sm_rootzone, :sm_rootzone_pctl, :sm_rootzone_wetness, :sm_surface, :sm_surface_wetness, 
    :snow_depth, :snow_mass, :snow_melt_flux, :snowfall_surface_flux, :soil_temp_layer1, 
    :soil_temp_layer2, :soil_temp_layer3, :soil_temp_layer4, :soil_temp_layer5, :soil_temp_layer6, 
    :soil_water_infiltration_flux, :specific_humidity_lowatmmodlay, :surface_pressure, :surface_temp, 
    :temp_lowatmmodlay, :vegetation_greenness_fraction, :windspeed_lowatmmodlay
)

if isfile(path1) && isfile(path2)
    # We need to wrap HDF5 for SMAP, as h5 file may not be SMAP files
    @testset "SMAPhdf5 wrapper" begin
        HDF5.h5open(path1) do f
            ds= GeoData.SMAPhdf5(f)
            @test keys(ds) == layerkeys(ds) == smapkeys
            @test dims(ds) isa Tuple{<:X,<:Y}
        end
    end

    @testset "geoarray" begin
        smaparray = geoarray(path1)

        @testset "open" begin
            @test all(open(A -> A[Y=1], smaparray) .=== smaparray[:, 1])
        end

        @testset "read" begin
            A = read(smaparray)
            @test A isa GeoArray
            @test parent(A) isa Array
        end

        @testset "array properties" begin
            @test smaparray  isa GeoArray
        end

        @testset "dimensions" begin
            @test ndims(smaparray) == 2
            HDF5.h5open(path1) do ds
                @test size(smaparray) == length.(dims(smaparray)) == size(ds["Geophysical_Data/baseflow_flux"])
            end
            @test dims(smaparray) isa Tuple{<:X,<:Y}
            @test refdims(smaparray) == ()
            # TODO detect the time span, and make it Regular
            modes = (
                Mapped(Ordered(), Irregular((-180.0f0, 180.0f0)), Intervals(Center()), GeoData.SMAPCRS, EPSG(4326)),
                Mapped(Ordered(ReverseIndex(), ReverseArray(), ForwardRelation()), Irregular((-85.04456f0, 85.04456f0)), Intervals(Center()), GeoData.SMAPCRS, EPSG(4326)),
            )
            @test val.(span(smaparray)) == val.(span.(modes))
            @test typeof(mode(smaparray)) == typeof(modes)
            @test bounds(smaparray) == ((-180.0f0, 180.0f0), (-85.04456f0, 85.04456f0))
        end

        @testset "other fields" begin @test missingval(smaparray) == -9999.0
            @test metadata(smaparray) isa Metadata{SMAPfile}
            @test name(smaparray) == :baseflow_flux
        end

        @testset "indexing" begin
            @test smaparray[Ti(1)] isa GeoArray{<:Any,2}
            @test smaparray[Y(1), Ti(1)] isa GeoArray{<:Any,1}
            @test smaparray[X(1), Ti(1)] isa GeoArray{<:Any,1}
            @test smaparray[X(1), Y(1)] == -9999.0
            @test smaparray[X(30), Y(30)] isa Float32
            # Russia
            @test smaparray[X(50), Y(100)] == -9999.0 
            # Alaska
            @test smaparray[Y(Near(64.2008)), X(Near(149.4937))] == 0.0f0
        end

        @testset "selectors" begin
            a = smaparray[X(Near(21.0)), Y(Between(50, 52))]
            index(smaparray, Y)
            indexorder(smaparray, Y)
            @test bounds(a) == ((50.08451f0, 51.977905f0),)
            x = smaparray[Lon(Near(150)), Lat(Near(30))]
            @test x isa Float32
            dimz = Lon(Between(-180.0, 180)), Lat(Between(-90, 90)) 
            @test size(smaparray[dimz...]) == (3856, 1624)
            @test index(smaparray[dimz...]) == index(smaparray)
            nca = smaparray[Lat(Between(-80, -25)), Lon(Between(0, 180))]
        end

        @testset "conversion to memory back GeoArray" begin
            geoA = read(smaparray[X(1:50), Y(20:20)])
            @test size(geoA) == (50, 1)
            @test eltype(geoA) <: Union{Missing,Float32}
            @time geoA isa GeoArray{Float32,1}
            @test dims(geoA) isa Tuple{<:X,<:Y}
            @test metadata(geoA) == metadata(smaparray)
            @test missingval(geoA) == missingval(smaparray)
            @test name(geoA) == :baseflow_flux
        end

        @testset "save" begin
            @testset "to grd" begin
                # TODO save and load subset
                geoA = read(smaparray)
                @test size(geoA) == size(smaparray)
                filename = tempname() * ".grd"
                write(filename, geoA)
                saved = read(geoarray(filename)[Band(1)])
                @test size(saved) == size(geoA)
                @test missingval(saved) === missingval(geoA)
                @test map(metadata.(dims(saved)), metadata.(dims(geoarray))) do s, g
                    all(s .== g)
                end |> all
                @test GeoData.name(saved) == GeoData.name(geoA)
                @test all(mode.(dims(saved)) .!= mode.(dims(geoA)))
                # @test all(order.(dims(saved)) .== order.(dims(geoA)))
                # @test all(typeof.(span.(dims(saved))) .== typeof.(span.(dims(geoA))))
                # @test all(val.(span.(dims(saved))) .== val.(span.(dims(geoA))))
                # @test all(sampling.(dims(saved)) .== sampling.(dims(geoA)))
                # @test typeof(dims(saved)) <: typeof(dims(geoA))
                # @test val(dims(saved)[2]) == val(dims(geoA)[2])
                # @test all(val.(dims(saved)) .== val.(dims(geoA)))
                # @test all(data(saved) .=== data(geoA))
                # @test saved isa typeof(geoA)
                # TODO test crs
            end
            @testset "to gdal" begin
                gdalfilename = tempname() * ".tif"
                write(gdalfilename, smaparray)
                gdalarray = geoarray(gdalfilename; mappedcrs=EPSG(4326))
                # These come out with slightly different format
                # @test rs(gdalarray) == crs(smaparray)
                @test_broken mappedbounds(gdalarray) == (bounds(smaparray)..., (1, 1))
                # Tiff locus = Start, Netcdf locus = Center
                @test Float32.(mappedindex(gdalarray, Y)) ≈ index(smaparray, Y)
                @test mappedindex(gdalarray, X) ≈ mappedindex(smaparray, X)
                @test all(gdalarray .== smaparray)
            end
            @testset "to netcdf" begin
                ncdfilename = tempname() * ".nc"
                write(ncdfilename, smaparray)
                saved = geoarray(ncdfilename)
                @test_broken bounds(saved) == bounds(smaparray)
                @test index(saved, Y) == reverse(index(smaparray, Y))
                @test index(saved, X) == index(smaparray, X)
            end
        end

        @testset "show" begin
            sh = sprint(show, MIME("text/plain"), smaparray)
            @test occursin("GeoArray", sh)
            @test occursin("Y", sh)
            @test occursin("X", sh)
        end

        @testset "plot" begin
            smaparray[Ti(1:3:12)] |> plot
            smaparray[Ti(1)] |> plot
            smaparray[Y(100), Ti(1)] |> plot
        end

    end

    @testset "stack" begin
        smapstack = stack(path1)

        @testset "read" begin
            st = read(smapstack)
            @test st isa GeoStack
            @test st.data isa NamedTuple
            @test first(st.data) isa Array
        end

        @testset "conversion to GeoArray" begin
            smaparray = smapstack[:soil_temp_layer1];
            @test smaparray isa GeoArray{Float32,2}
            @test dims(smaparray) isa Tuple{<:X{<:Array{Float32,1}}, <:Y{<:Array{Float32,1}}}
            @test span(smaparray) isa Tuple{Irregular{Tuple{Float32,Float32}},Irregular{Tuple{Float32,Float32}}}
            @test span(smaparray) == (Irregular((-180.0f0, 180.0f0)), Irregular((-85.04456f0, 85.04456f0)))
            @test bounds(smaparray) == ((-180.0f0, 180.0f0), (-85.04456f0, 85.04456f0))
            @test index(smaparray) isa Tuple{Vector{Float32},Vector{Float32}}
            @test refdims(smaparray) isa Tuple{<:Ti}
            @test missingval(smaparray) == -9999.0
            @test smaparray[1, 1] == -9999.0
            @test name(smaparray) == :soil_temp_layer1
            dt = DateTime(2016, 1, 1, 22, 30)
            step_ = Hour(3)
            @test refdims(smapstack) ==
                (Ti(dt:step_:dt; mode=Sampled(Ordered(), Regular(step_), Intervals(Start()))),)
            # Currently empty
            @test metadata(smaparray) isa Metadata{SMAPfile}
        end

        @testset "conversion to GeoStack" begin
            # Stack Constructors
            # This uses too much ram! There is a lingering memory leak in HDF5.
            # geostack = GeoStack(stack)
            # @test Symbol.(Tuple(keys(stack))) == keys(geostack)
            geostack = GeoStack(smapstack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
            @test keys(geostack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
        end

        if VERSION > v"1.1-"
            @testset "copy" begin
                geoA = zero(smapstack[:soil_temp_layer1])
                @test geoA isa GeoArray
                @test geoA != smapstack[:soil_temp_layer1]
                copy!(geoA, smapstack, :soil_temp_layer1)
                @test geoA == smapstack[:soil_temp_layer1]
            end
        end

        @testset "window" begin
            windowedstack = stack(path1; window=(Y(1:5), X(1:5)))
            @test window(windowedstack) == (Y(1:5), X(1:5))
            windowedarray = windowedstack[:soil_temp_layer1];
            @test size(windowedarray) == (5, 5)
            @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2] == -9999.0
            windowedstack = stack(path1; window=(Y(1:5), X(1:5), Ti(1:1)))
            windowedarray = windowedstack[:soil_temp_layer1];
            @test size(windowedarray) == (5, 5)
            @test windowedarray[1:3, 2:2, 1] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2, 1] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2, 1] == -9999.0
            windowedstack = stack(path1; window=(Ti(1),))
            windowedarray = windowedstack[:soil_temp_layer1];
            @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2] == -9999.0
        end

        @testset "show" begin
            sh1 = sprint(show, MIME("text/plain"), smapstack[:soil_temp_layer1])
            @test occursin("GeoArray", sh1)
            @test occursin("Y", sh1)
            @test occursin("X", sh1)
            @test occursin("Ti", sh1)
        end

    end

    @testset "series" begin
        ser = smapseries([path1, path2])
        val.(dims(ser))
        @test ser[1] isa GeoStack
        @test first(bounds(ser, Ti)) == DateTime(2016, 1, 1, 22, 30)
        @test last(bounds(ser, Ti)) == DateTime(2016, 1, 3, 1, 30)
        modified_series = modify(Array, ser)
        stackkeys = keys(modified_series[1])
        @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{stackkeys,<:Tuple{<:Array{Float32,2,},Vararg}}}}

        @testset "read" begin
            geoseries = read(ser)
            @test geoseries isa GeoSeries{<:GeoStack}
            @test geoseries.data isa Vector{<:GeoStack}
            @test first(geoseries.data[1].data) isa Array 
        end

        @testset "show" begin
            sh = sprint(show, MIME("text/plain"), ser)
            # Test but don't lock this down too much
            @test occursin("GeoSeries", sh)
            @test occursin("Ti", sh)
        end
    end
end
