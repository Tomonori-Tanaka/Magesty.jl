"""
Unit tests for the `Magesty.IncarIO` INCAR text reader/writer. Covers value
typing on parse, MAGMOM/M_CONSTR vector formatting, and a parse -> write ->
parse round trip that must reproduce the parameter dictionary exactly.
"""

using Magesty
using Test
using DataStructures: OrderedDict

using Magesty.IncarIO: parse_incar, write_incar

const INCAR_FIXTURE = joinpath(@__DIR__, "fixtures", "incar", "INCAR")

@testset "parse_incar typing" begin
    d = parse_incar(INCAR_FIXTURE)

    @test d[:SYSTEM] == "test system"
    @test d[:ENCUT] === 520.0          # trailing ".0" -> Float64
    @test d[:ISMEAR] === 0             # integral -> Int
    @test d[:SIGMA] === 0.05
    @test d[:LSORBIT] === true         # .TRUE. -> Bool

    # MAGMOM / M_CONSTR are always Vector{Float64}, length a multiple of 3.
    @test d[:MAGMOM] isa Vector{Float64}
    @test d[:M_CONSTR] isa Vector{Float64}
    @test length(d[:MAGMOM]) == 12
    @test d[:MAGMOM] == [0.0, 0.0, 3.0, 3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, -3.0]
    @test d[:M_CONSTR] == d[:MAGMOM]

    # Insertion order is preserved.
    @test collect(keys(d)) == [:SYSTEM, :ENCUT, :ISMEAR, :SIGMA, :LSORBIT, :MAGMOM, :M_CONSTR]

    missing_path = joinpath(@__DIR__, "fixtures", "incar", "DOES_NOT_EXIST")
    @test_throws ArgumentError parse_incar(missing_path)
end

@testset "n*v repeat-pattern expansion" begin
    # "n*v" expands to n copies of v (VASP MAGMOM zero-padding idiom).
    tmp = tempname()
    try
        write(tmp, "MAGMOM = 3*0.0 0.0 0.0 2.0\n")
        d = parse_incar(tmp)
        @test d[:MAGMOM] == [0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    finally
        isfile(tmp) && rm(tmp)
    end
end

@testset "parse -> write -> parse round trip" begin
    d = parse_incar(INCAR_FIXTURE)
    for wrap in (false, true)
        tmp = tempname()
        try
            write_incar(tmp, d; wrap_vectors = wrap)
            d2 = parse_incar(tmp)
            @test collect(keys(d2)) == collect(keys(d))
            for k in keys(d)
                @test d2[k] == d[k]
            end
        finally
            isfile(tmp) && rm(tmp)
        end
    end
end

@testset "MAGMOM/M_CONSTR formatting" begin
    d = OrderedDict{Symbol, Any}(:MAGMOM => [0.0, 0.0, 1.5, 0.0, 0.0, -1.5])
    tmp = tempname()
    try
        write_incar(tmp, d)
        text = read(tmp, String)
        # 3 components per atom, double space between atoms, 9 decimals.
        @test occursin("MAGMOM = 0.000000000 0.000000000 1.500000000  " *
                       "0.000000000 0.000000000 -1.500000000", text)
    finally
        isfile(tmp) && rm(tmp)
    end
end

println("All IncarIO tests passed.")
