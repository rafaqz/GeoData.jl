using GeoData, RasterDataSources, Test, Dates, NCDatasets, ArchGDAL

# Too big to test on CI
if !haskey(ENV, "CI")
    @testset "load WorldClim Weather" begin
        # Weather time-series
        dates = (Date(2001), Date(2002))
        ser = series(WorldClim{Weather}, (:prec,); date=dates, window=(Lat(600:1900), Lon(600:1900)))
        ser[Date(2001, 1)][:prec]
        # Select Australia, using regular lat/lon selectors
        A = geoarray(WorldClim{Weather}, :prec; date=DateTime(2001, 05), mappedcrs=EPSG(4326))
    end
end

@testset "load WorldClim Climate" begin
    # Weather time-series
    ser = series(WorldClim{Climate}, (:prec,); res="10m", month=Jan:March);
    ser[Jan][:prec] 
    # Select Australia, using regular lat/lon selectors
    A = geoarray(WorldClim{Climate}, :prec; month=1, mappedcrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))]
    st = stack(WorldClim{BioClim}, (1, 2))
    st[:BIO1]
    @test st isa GeoStack
    @test A isa GeoArray
end

@testset "load WorldClim BioClim" begin
    A = geoarray(WorldClim{BioClim}, 1; mappedcrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))]
    @test A isa GeoArray
end

@testset "load CHELSA BioClim" begin
    A = geoarray(CHELSA{BioClim}, 1; mappedcrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))]
    st = stack(CHELSA{BioClim}, (1, 2))
    @test A isa GeoArray
    @test st isa GeoStack
    @test st[:BIO2] isa GeoArray
end

@testset "load EarthEnv HabitatHeterogeneity" begin
    A = geoarray(EarthEnv{HabitatHeterogeneity}, :cv; mappedcrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))] 
    st = stack(EarthEnv{HabitatHeterogeneity}, (:cv, :evenness))
    @test A isa GeoArray
    @test st isa GeoStack
    @test st[:evenness] isa GeoArray
end

@testset "load EarthEnv LandCover" begin
    A = geoarray(EarthEnv{LandCover}, 2; mappedcrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))]
    @test A isa GeoArray
end

@testset "load ALWB" begin
    A = geoarray(ALWB{Deciles,Day}; date=DateTime(2019, 10, 19))
    @test crs(A) == EPSG(4326)
    A = geoarray(ALWB{Values,Day}, :ss_pct; date=DateTime(2019, 10, 19))
    @test crs(A) == EPSG(4326)
    st = stack(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=DateTime(2019, 10, 19))
    @test crs(st[:s0_pct]) == EPSG(4326)
    dates = DateTime(2019, 10, 19):Day(1):DateTime(2019, 11, 20)
    s = series(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=dates)
    @test A isa GeoArray
    @test st isa GeoStack
    @test s isa GeoSeries
end

# Obscure .Z format may not work on windows
if Sys.islinux()
    @testset "load AWAP" begin
        A = geoarray(AWAP, :rainfall; date=DateTime(2019, 10, 19))
        st = stack(AWAP; date=DateTime(2019, 10, 19), resize=crop)
        dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
        s = series(AWAP; date=dates, resize=crop);
        s[DateTime(2019, 10, 1)][:solar]
        @test A isa GeoArray
        @test st isa GeoStack
        @test s isa GeoSeries
    end
end
