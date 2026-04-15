"""
Unit tests for ExtXYZWriter.

Tests are self-contained — no fixture files required.
Run from repo root:
    julia tools/test/test_ExtXYZ.jl
"""

include(joinpath(@__DIR__, "../ExtXYZ.jl"))
using .ExtXYZWriter
using Test
using LinearAlgebra

# ── helpers ───────────────────────────────────────────────────────────────────

"""Parse the extxyz header key=value pairs into a Dict{String,String}."""
function parse_header(line::AbstractString)::Dict{String, String}
    d = Dict{String, String}()
    # Tokenise: quoted values or bare tokens
    i = 1
    while i <= length(line)
        m = match(r"(\w+)=(\"[^\"]*\"|\S+)", line[i:end])
        m === nothing && break
        key = m.captures[1]
        val = strip(m.captures[2], '"')
        d[key] = val
        i += i + m.offset + length(m.match) - 2
    end
    return d
end

# ── minimal frame ─────────────────────────────────────────────────────────────

@testset "minimal frame (no forces/stress/extras)" begin
    lat = Matrix{Float64}(I, 3, 3) .* 5.0
    frame = AtomFrame(
        num_atoms   = 2,
        lattice     = lat,
        pbc         = [true, true, false],
        species     = ["Fe", "Fe"],
        positions   = [0.0 2.5; 0.0 2.5; 0.0 2.5],
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    lines = split(String(take!(buf)), "\n"; keepempty=false)

    @test lines[1] == "2"                                 # atom count
    hdr = lines[2]
    @test occursin("Lattice=", hdr)
    @test occursin("Properties=species:S:1:pos:R:3", hdr)
    @test occursin("pbc=\"T T F\"", hdr)
    @test !occursin("forces", hdr)
    @test !occursin("energy_free", hdr)
    @test !occursin("stress", hdr)
    @test length(lines) == 4                              # header + 2 atoms
end

# ── decimal alignment (space flag) ────────────────────────────────────────────

@testset "space-flag decimal alignment" begin
    lat = Matrix{Float64}(I, 3, 3) .* 3.0
    frame = AtomFrame(
        num_atoms   = 2,
        lattice     = lat,
        species     = ["Fe", "Fe"],
        positions   = [1.0 -1.0; 0.0 0.0; 0.0 0.0],
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    lines = split(String(take!(buf)), "\n"; keepempty=false)

    # The space flag ("%") inserts a leading space for positive numbers and "-"
    # for negative, so consecutive fields are separated by two spaces (space +
    # delimiter) for positive values and one space + "-" for negative.
    pos_line = lines[3]   # first atom: x=1.0
    neg_line = lines[4]   # second atom: x=-1.0
    # After "Fe " the x-position field starts: "  1." for positive (space flag + delimiter space)
    @test occursin("  1.", pos_line)
    # For negative: " -1."
    @test occursin(" -1.", neg_line)
end

# ── scalar extra_per_atom (ncols=1) ───────────────────────────────────────────

@testset "scalar extra_per_atom → Properties ncols=1" begin
    lat = Matrix{Float64}(I, 3, 3) .* 4.0
    rwigs = reshape([1.3, 1.5], 1, 2)   # 1×2 matrix
    frame = AtomFrame(
        num_atoms      = 2,
        lattice        = lat,
        species        = ["Mn", "Ir"],
        positions      = zeros(3, 2),
        extra_per_atom = ["rwigs" => rwigs],
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    lines = split(String(take!(buf)), "\n"; keepempty=false)
    hdr = lines[2]

    @test occursin("rwigs:R:1", hdr)
    # Each atom line should have: species + 3 pos + 1 rwigs = 5 tokens
    @test length(split(lines[3])) == 5
    @test length(split(lines[4])) == 5
end

# ── vector extra_per_atom (ncols=3) ───────────────────────────────────────────

@testset "vector extra_per_atom → Properties ncols=3" begin
    lat = Matrix{Float64}(I, 3, 3) .* 4.0
    magmom = [1.0 -1.0; 0.0 0.0; 0.5 -0.5]   # 3×2
    frame = AtomFrame(
        num_atoms      = 2,
        lattice        = lat,
        species        = ["Fe", "Fe"],
        positions      = zeros(3, 2),
        forces         = zeros(3, 2),
        extra_per_atom = ["MAGMOM_smoothed" => magmom],
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    lines = split(String(take!(buf)), "\n"; keepempty=false)
    hdr = lines[2]

    @test occursin("MAGMOM_smoothed:R:3", hdr)
    # Each atom line: species + 3 pos + 3 forces + 3 magmom = 10 tokens
    @test length(split(lines[3])) == 10
end

# ── energy and stress in header ───────────────────────────────────────────────

@testset "energy_free, energy_zero, stress in header" begin
    lat = Matrix{Float64}(I, 3, 3) .* 3.0
    stress = [1.0 0.1 0.2; 0.1 2.0 0.3; 0.2 0.3 3.0]  # 3×3 eV/Å³
    frame = AtomFrame(
        num_atoms    = 1,
        lattice      = lat,
        species      = ["Fe"],
        positions    = zeros(3, 1),
        energy_free  = -10.123,
        energy_zero  = -10.456,
        stress       = stress,
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    hdr = split(String(take!(buf)), "\n")[2]

    @test occursin("energy_free=", hdr)
    @test occursin("energy_zero=", hdr)
    # Stress Voigt: xx yy zz yz xz xy
    @test occursin("stress=", hdr)
    # Parse the stress values from the header
    m = match(r"stress=\"([^\"]+)\"", hdr)
    @test m !== nothing
    voigt = parse.(Float64, split(m.captures[1]))
    @test length(voigt) == 6
    @test voigt[1] ≈ 1.0    # xx
    @test voigt[2] ≈ 2.0    # yy
    @test voigt[3] ≈ 3.0    # zz
    @test voigt[4] ≈ 0.3    # yz
    @test voigt[5] ≈ 0.2    # xz
    @test voigt[6] ≈ 0.1    # xy
end

# ── comment field ─────────────────────────────────────────────────────────────

@testset "comment field in header" begin
    lat = Matrix{Float64}(I, 3, 3) .* 3.0
    frame = AtomFrame(
        num_atoms = 1,
        lattice   = lat,
        species   = ["Fe"],
        positions = zeros(3, 1),
        comment   = "VASP6.6.0, ENCUT=300",
    )

    buf = IOBuffer()
    write_extxyz(buf, frame)
    hdr = split(String(take!(buf)), "\n")[2]
    @test occursin("comment=\"VASP6.6.0, ENCUT=300\"", hdr)
end

# ── file output ───────────────────────────────────────────────────────────────

@testset "write to file and read back" begin
    lat = Matrix{Float64}(I, 3, 3) .* 4.0
    frame = AtomFrame(
        num_atoms = 2,
        lattice   = lat,
        species   = ["Fe", "Rh"],
        positions = [0.0 2.0; 0.0 2.0; 0.0 2.0],
    )

    tmp = tempname() * ".extxyz"
    try
        write_extxyz(tmp, frame)
        lines = readlines(tmp)
        @test lines[1] == "2"
        @test occursin("Lattice=", lines[2])
        @test lines[3][1:2] == "Fe"
        @test lines[4][1:2] == "Rh"
    finally
        isfile(tmp) && rm(tmp)
    end
end

println("All ExtXYZWriter tests passed.")
