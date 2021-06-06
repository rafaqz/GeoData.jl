using GeoData, DimensionalData, Test, Statistics, Dates
using GeoData: data, window

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X((10, 100)), Y((-50, 50))
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
nme = :test
mval = -9999.0
meta = NoMetadata()

# Formatting only occurs in shorthand constructors
ga1 = GeoArray(data1, dims1; refdims=refdimz, name=nme, metadata=meta, missingval=mval)
ga2 = GeoArray(data2, dims2)

st = GeoStack(ga1, ga2; keys=(:ga1, :ga2))

@testset "stack layers" begin
    @test length(st) == 2
    @test first(st) == ga1
    @test last(st) == ga2
    @test DimensionalData.layers(st) isa NamedTuple
    @test st[:ga1] == ga1
    @test st[:ga2] == ga2
    @test data(st[:ga1]) == data1
    @test data(st[:ga1]) isa Array{Float64,2}
    @test keys(st) == (:ga1, :ga2)
    @test haskey(st, :ga1)
    @test names(st) == (:ga1, :ga2)
    @test collect(values(st)) == [ga1, ga2]
    # We can iterate/splat stacks as GeArrays
    @test st === GeoStack(st...)
end

@testset "st fields " begin
    @test DimensionalData.layerdims(st, :ga1) == DimensionalData.formatdims(data1, dims1)
    @test window(st) == nothing
    @test metadata(st) == NoMetadata()
    @test metadata(st, :ga1) == NoMetadata()
end

@testset "indexing" begin
    # Indexing the st is the same as indexing its child array
    a = st[:ga1][X(2:4), Y(5:6)]
    @test a == st[:ga1, X(2:4), Y(5:6)]

    @inferred st[:ga1][X=2:4, Y=5:6]
    # FIXME: This isn't inferred, the constants don't propagate like they
    # do in the above call. Probably due to the anonymous wrapper function.
    @test_broken @inferred st[:ga1, X(2:4), Y(5:6)] isa GeoArray

    # Getindex for a whole st of new GeoArrays
    a = st[X=2:4, Y=5:6]
    @test a isa GeoStack
    @test a[:ga1] isa GeoArray
    @test data(a[:ga1]) isa Array
    @test a[:ga1] == data1[2:4, 5:6]
    @test a[:ga2] == data2[2:4, 5:6, 1:1]

    @testset "select new arrays for the whole st" begin
        s = st[Y=Between(-10, 10.0), Ti=At(DateTime(2019))]
        @test s isa GeoStack
        @test s[:ga1] isa GeoArray
        @test data(s[:ga1]) isa Array
        @test s[:ga1] == data1[:, 5:7]
        @test s[:ga2] == data2[:, 5:7, 1]
        @test dims(s[:ga2]) == (X(LinRange(10.0, 100.0, 10); mode=Sampled(Ordered(), Regular(10.0), Points())),
                                Y(LinRange(-10.0, 10.0, 3); mode=Sampled(Ordered(), Regular(10.0), Points())))
        @test refdims(s[:ga2]) == (Ti(DateTime(2019); mode=Sampled(Ordered(), Irregular(), Points())),)
        @test ismissing(missingval(s, :ga2)) && ismissing(missingval(s[:ga2]))
    end

    @testset "select views of arrays for the whole st" begin
        sv = view(st, Y=Between(-4.0, 27.0), Ti=At(DateTime(2019)))
        @test sv isa GeoStack
        @test sv[:ga1] isa GeoArray
        @test data(sv[:ga1]) isa SubArray
        @test sv[:ga1] == data1[:, 6:8]
        @test sv[:ga2] == data2[:, 6:8, 1]
        @test dims(sv[:ga2]) == (X(LinRange(10.0, 100.0, 10); mode=Sampled(Ordered(), Regular(10.0), Points())),
                                 Y(LinRange(0.0, 20.0, 3); mode=Sampled(Ordered(), Regular(10.0), Points())))
        @test refdims(sv[:ga2])[1] == (Ti(DateTime(2019); mode=Sampled(Ordered(), Irregular(), Points())),)[1]
        # Stack of view-based GeoArrays
        v = view(st, X(2:4), Y(5:6))
        @inferred view(st, X(2:4), Y(5:6))
        @test v isa GeoStack
        @test v[:ga1] isa GeoArray
        @test data(v[:ga1]) isa SubArray
        @test v[:ga1] == view(data1, 2:4, 5:6)
        @test v[:ga2] == view(data2, 2:4, 5:6, 1:1)
    end

end

@testset "subset st with specific key(s)" begin
    s1 = GeoStack(st; keys=(:ga2,))
    @test keys(s1) == (:ga2,)
    @test length(values(s1)) == 1
    s2 = GeoStack(st; keys=(:ga1, :ga2))
    @test keys(s2) == (:ga1, :ga2)
    @test length(values(s2)) == 2
end

@testset "concatenate stacks" begin
    dims1b = X((110, 200)), Y((-50, 50))
    dims2b = (dims1b..., Ti([DateTime(2019)]))
    stack_a = GeoStack((ga1=ga1, ga2=ga2))
    stack_b = GeoStack((ga1=GeoArray(data1 .+ 10, dims1b), ga2=GeoArray(data2 .+ 20, dims2b)))
    catstack = cat(stack_a, stack_b; dims=X())
    @test size(first(catstack)) == (20, 11)
    @test val(dims(catstack, X)) ≈ 10.0:10.0:200.0
    #@test step(dims(first(catstack), X())) == 10.0
    @test DimensionalData.bounds(dims(first(catstack), X)) == (10.0, 200.0)
    @test catstack[:ga1][Y(1)] == 1.0:20.0
    @test catstack[:ga2][Y(1), Ti(1)] == 2.0:2.0:40.0
    catstack = cat(stack_a, stack_b; dims=(X(), Y()))
    @test size(first(catstack)) == (20, 22)
end

@testset "show" begin
    sh = sprint(show, st)
    # Test but don't lock this down too much
    @test occursin("GeoStack", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end
