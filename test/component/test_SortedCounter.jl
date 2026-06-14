using .SortedCounters

@testset "SortedCounter" begin
    # Empty construction
    sc = SortedCounter{Int}()
    @test isempty(sc)
    @test length(sc) == 0
    @test collect(sc) == Int[]

    # push! single
    push!(sc, 3)
    push!(sc, 1)
    push!(sc, 3)
    @test length(sc) == 2
    @test !isempty(sc)
    @test sc.counts[1] == 1
    @test sc.counts[3] == 2
    @test collect(sc) == [1, 3]
    @test sc[1] == 1
    @test sc[2] == 3

    # push! with count
    push!(sc, 2, 5)
    @test length(sc) == 3
    @test sc.counts[2] == 5
    @test collect(sc) == [1, 2, 3]

    # push! with count on existing key
    push!(sc, 3, 10)
    @test sc.counts[3] == 12

    # == compares counts
    sc2 = SortedCounter{Int}()
    push!(sc2, 2, 5)
    push!(sc2, 1)
    push!(sc2, 3, 12)
    @test sc == sc2

    # isless compares sorted keys
    sc_a = SortedCounter{Int}()
    push!(sc_a, 5)
    sc_b = SortedCounter{Int}()
    push!(sc_b, 1)
    @test sc_b < sc_a

    # copy independence
    sc_orig = SortedCounter{Int}()
    push!(sc_orig, 7)
    collect(sc_orig)  # force cache refresh
    sc_copy = copy(sc_orig)
    push!(sc_copy, 9)
    @test !haskey(sc_orig.counts, 9)
    @test haskey(sc_copy.counts, 9)
    @test length(sc_orig) == 1
    @test length(sc_copy) == 2

    # Vector{Int} keys (the irreducible cluster use case)
    scv = SortedCounter{Vector{Int}}()
    push!(scv, [1, 2, 3])
    push!(scv, [1, 2, 4])
    push!(scv, [1, 2, 3])
    @test length(scv) == 2
    @test scv.counts[[1, 2, 3]] == 2
    @test scv.counts[[1, 2, 4]] == 1
    @test collect(scv) == [[1, 2, 3], [1, 2, 4]]

    # Cache invalidation across alternating push! and iterate
    sc_iv = SortedCounter{Int}()
    push!(sc_iv, 5)
    @test collect(sc_iv) == [5]   # first refresh
    push!(sc_iv, 2)               # invalidates
    @test collect(sc_iv) == [2, 5]
    push!(sc_iv, 5)               # increment existing — no invalidation but value changes
    @test collect(sc_iv) == [2, 5]
    @test sc_iv.counts[5] == 2

    # .counts is exposed for O(1) lookup with default
    @test get(sc_iv.counts, 999, 0) == 0
    @test get(sc_iv.counts, 5, 0) == 2
end
