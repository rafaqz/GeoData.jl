using GeoData, Test, Dates

# GeoSeries from GeoArray/GeoStack components

data1 = [1 2 3 4
         5 6 7 8]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = X([30, 40]), Y((-10, 20))
ga1 = GeoArray(data1, dimz)
ga2 = GeoArray(data2, dimz)
ga1a = GeoArray(data3, dimz)
ga2a = GeoArray(data4, dimz)
stack1 = GeoStack(ga1, ga2; keys=(:ga1, :ga2))
stack2 = GeoStack(ga1a, ga2a; keys=(:ga1, :ga2))
dates =[DateTime(2017), DateTime(2018)]
ser = GeoSeries([stack1, stack2], (Ti(dates),))
@test issorted(dates)

@testset "getindex returns the currect types" begin
    @test ser[Ti(1)] isa GeoStack{<:NamedTuple}
    @test ser[Ti(1)][:ga2] isa GeoArray{Int,2}
    @test ser[Ti(1)][:ga2, 1, 1] isa Int
    @test ser[Ti(1)][:ga2][1, 1] isa Int
end

@testset "properties" begin
    @test refdims(ser) === ()
    # Should these be real fields? what is the use-case?  @test metadata(ser) === nothing @test name(ser) === ""
    @test label(ser) === ""
end

@testset "getindex returns the currect results" begin
    @test ser[Ti(Near(DateTime(2017)))][:ga1][X(1), Y(3)] === 3
    @test ser[Ti(At(DateTime(2017)))][:ga1, X(1), Y(3)] === 3
    @test ser[Ti(At(DateTime(2018)))][:ga2][X(2), Y(4)] === 32
    @test ser[Ti(At(DateTime(2018)))][:ga2, X(2), Y(4)] === 32
    @test ser[Ti(1)][:ga1, X(1), Y(2)] == 2
    @test ser[Ti(1)][:ga2, X(2), Y(3:4)] == [14, 16] 
end

@testset "getindex is type stable all the way down" begin
    # @inferred ser[Ti(At(DateTime(2017)))][:ga1, X(1), Y(2)]
    @inferred ser[Ti(1)][:ga1][X(1), Y(2)]
    # @inferred ser[Ti(1)][:ga1, X(1), Y(2:4)]
    @inferred ser[Ti(1)][:ga1][X(1), Y(2:4)]
    # @inferred ser[1][:ga1, X(1:2), Y(:)]
    @inferred ser[1][:ga1][X(1:2), Y(:)]
end

@testset "lazy view windows" begin
    dimz = (Ti([DateTime(2017), DateTime(2018)]),)
    dat = [stack1, stack2]
    window_ = X(1:2), Y(3:4)
    ser = GeoSeries(dat, dimz; window=window_)
    st = ser[1]
    @test st[:ga1] == [3 4; 7 8]
    @test st[:ga1, 1, 2] == 4
end

@testset "setindex!" begin
    ser[1] = ser[2]
    @test typeof(ser[1]) == typeof(ser[2]) == eltype(parent(ser))
    typeof(parent(ser)[1]) == eltype(ser)
    ser[Ti(1)] = ser[Ti(2)]
    @test ser[Ti(1)] == ser[Ti(2)]
end

