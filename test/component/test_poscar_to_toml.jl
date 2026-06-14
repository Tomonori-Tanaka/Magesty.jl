"""
Regression tests for the `poscar_to_toml` API. The output must be
byte-identical to the golden TOML fixtures. The Direct-coordinate golden
file was generated from the POSCAR fixture; the Cartesian golden file was
derived analytically (it exercises the scaling-factor handling that the
Cartesian-to-direct conversion depends on). The `magesty vasp toml` CLI
subcommand is covered separately by the `MagestyCLI` package tests.
"""

using Magesty
using Test

@testset "poscar_to_toml: $(case)" for (case, poscar, golden) in (
        ("direct",    "POSCAR",           "expected.toml"),
        ("cartesian", "POSCAR_cartesian", "expected_cartesian.toml"),
    )
    dir         = joinpath(@__DIR__, "fixtures", "poscar")
    poscar_path = joinpath(dir, poscar)
    golden_text = read(joinpath(dir, golden), String)

    # API function: the returned text is byte-identical to the golden file.
    @test poscar_to_toml(poscar_path) == golden_text

    # The `output` keyword writes the same bytes to disk.
    api_outfile = tempname() * ".toml"
    try
        poscar_to_toml(poscar_path; output = api_outfile)
        @test read(api_outfile, String) == golden_text
    finally
        isfile(api_outfile) && rm(api_outfile)
    end
end

println("All poscar_to_toml tests passed.")
