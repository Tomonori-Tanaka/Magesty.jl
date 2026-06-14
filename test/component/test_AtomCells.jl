using Test
using Magesty

# `AtomCell` is an internal type: the `AtomCells` submodule exports it,
# but `Magesty` itself does not re-export, so the test reaches it
# through the canonical module path. Inside the package the type is
# used as an element of `Vector{AtomCell}` and as a `Set` / `Dict`
# key in `src/Clusters.jl`; the tests below pin the four `Base`
# methods (`==`, `isless`, `hash`, `show`) that back those usages.
const AtomCell = Magesty.AtomCells.AtomCell

@testset "AtomCells" begin
    # ---- equality (==) --------------------------------------------------

    @testset "== (equality)" begin
        a = AtomCell(1, 2)
        b = AtomCell(1, 2)

        # Reflexive and value-based: identity does not matter.
        @test a == a
        @test a == b
        @test b == a  # symmetric

        # Differing atom or cell breaks equality.
        @test a != AtomCell(2, 2)
        @test a != AtomCell(1, 3)
        @test a != AtomCell(2, 3)
    end

    # ---- isless: lexicographic over (atom, cell) ------------------------

    @testset "isless (total order)" begin
        # Primary key: atom. Secondary key: cell.
        @test isless(AtomCell(1, 5), AtomCell(2, 0))  # atom dominates
        @test isless(AtomCell(1, 1), AtomCell(1, 2))  # cell tiebreaks
        @test !isless(AtomCell(2, 0), AtomCell(1, 5)) # asymmetric
        @test !isless(AtomCell(1, 2), AtomCell(1, 1))
        @test !isless(AtomCell(1, 1), AtomCell(1, 1)) # strict

        # Antisymmetry on every distinct pair: exactly one of (a, b) /
        # (b, a) is true.
        pairs = [
            (AtomCell(1, 1), AtomCell(1, 2)),
            (AtomCell(1, 5), AtomCell(2, 0)),
            (AtomCell(2, 5), AtomCell(2, 6)),
        ]
        for (x, y) in pairs
            @test isless(x, y) != isless(y, x)
        end

        # `sort` produces the lexicographic order, which is exactly the
        # downstream `Clusters.jl` "sorted AtomCell order" requirement.
        input = [AtomCell(2, 0), AtomCell(1, 5), AtomCell(1, 1), AtomCell(2, 0)]
        expected = [AtomCell(1, 1), AtomCell(1, 5), AtomCell(2, 0), AtomCell(2, 0)]
        @test sort(input) == expected
    end

    # ---- hash: consistent with == --------------------------------------

    @testset "hash (consistent with ==)" begin
        a = AtomCell(3, 7)
        b = AtomCell(3, 7)

        # The `==`/`hash` contract: equal values must hash to the same
        # UInt under every seed.
        @test hash(a) == hash(b)
        @test hash(a, UInt(0)) == hash(b, UInt(0))
        @test hash(a, UInt(42)) == hash(b, UInt(42))

        # Hash is a pure function: repeated calls return the same UInt.
        @test hash(a) == hash(a)
        @test hash(a, UInt(0)) == hash(a, UInt(0))

        # The seed must mix in. If the seed were ignored, `AtomCell`
        # would not be safe as a nested `Dict` key.
        @test hash(a, UInt(0)) != hash(a, UInt(1))
    end

    # ---- show: documented format does not drift -------------------------

    @testset "show" begin
        # The docstring example reads `(atom: 1, cell: 2)`; pin it so
        # format-tweaking PRs see a clear test failure.
        @test repr(AtomCell(1, 2)) == "(atom: 1, cell: 2)"
        @test repr(AtomCell(0, -1)) == "(atom: 0, cell: -1)"
    end

    # ---- cross-method integration: Set / Dict --------------------------

    @testset "Set / Dict usage" begin
        # Set dedup uses (`hash`, `==`) together.
        s = Set([AtomCell(1, 2), AtomCell(1, 2), AtomCell(3, 4)])
        @test length(s) == 2
        @test AtomCell(1, 2) in s
        @test AtomCell(3, 4) in s
        @test !(AtomCell(2, 1) in s)

        # Dict round-trip: an equal-but-distinct AtomCell finds the
        # same value, and overwrites it in place rather than adding a
        # new entry.
        d = Dict{AtomCell, Int}()
        d[AtomCell(5, 6)] = 100
        @test d[AtomCell(5, 6)] == 100
        d[AtomCell(5, 6)] = 200
        @test length(d) == 1
        @test d[AtomCell(5, 6)] == 200
    end
end
