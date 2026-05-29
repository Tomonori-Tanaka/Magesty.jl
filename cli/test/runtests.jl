"""
Regression tests for the `magesty vasp` CLI subcommands (`extxyz`, `toml`,
`embset`). The command output must be byte-identical to the golden
fixtures, which live in the core package test tree
(`test/component/fixtures/`) and are reused here by resolving a path
relative to the repository root.
"""

using MagestyCLI
using Test

const FIXTURE_DIR = joinpath(@__DIR__, "..", "..", "test", "component", "fixtures")

@testset "magesty vasp extxyz: $(system)" for system in ("FeRh", "IrMn3")
	dir          = joinpath(FIXTURE_DIR, system)
	vasprun_path = joinpath(dir, "vasprun.xml")
	oszicar_path = joinpath(dir, "OSZICAR")
	golden_text  = read(joinpath(dir, "$(system).extxyz"), String)

	# CLI subcommand `magesty vasp extxyz`: exits 0 and writes the same bytes.
	cli_outfile = tempname() * ".extxyz"
	try
		exit_code = MagestyCLI.command_main(
			["vasp", "extxyz", vasprun_path, "--oszicar", oszicar_path, "--output", cli_outfile],
		)
		@test exit_code == 0
		@test read(cli_outfile, String) == golden_text
	finally
		isfile(cli_outfile) && rm(cli_outfile)
	end
end

@testset "magesty vasp toml" begin
	dir         = joinpath(FIXTURE_DIR, "poscar")
	poscar_path = joinpath(dir, "POSCAR")
	golden_text = read(joinpath(dir, "expected.toml"), String)

	# CLI subcommand `magesty vasp toml`: exits 0 and writes the same bytes.
	cli_outfile = tempname() * ".toml"
	try
		exit_code = MagestyCLI.command_main(
			["vasp", "toml", poscar_path, "--output", cli_outfile],
		)
		@test exit_code == 0
		@test read(cli_outfile, String) == golden_text
	finally
		isfile(cli_outfile) && rm(cli_outfile)
	end
end

@testset "magesty vasp embset" begin
	dir         = joinpath(FIXTURE_DIR, "outcar")
	outcar_path = joinpath(dir, "OUTCAR")
	golden_text = read(joinpath(dir, "expected.embset"), String)

	# CLI subcommand `magesty vasp embset`: exits 0 and writes the same bytes.
	cli_outfile = tempname()
	try
		exit_code = MagestyCLI.command_main(
			["vasp", "embset", outcar_path, "--output", cli_outfile],
		)
		@test exit_code == 0
		@test read(cli_outfile, String) == golden_text
	finally
		isfile(cli_outfile) && rm(cli_outfile)
	end
end

@testset "magesty sunny script" begin
	model_path = joinpath(@__DIR__, "..", "..", "test", "integration", "dimer", "dimer_dmi.xml")

	# CLI subcommand `magesty sunny script`: exits 0 and writes a parseable
	# Sunny.jl script.
	cli_outfile = tempname() * ".jl"
	try
		exit_code = MagestyCLI.command_main(
			["sunny", "script", model_path, "--output", cli_outfile],
		)
		@test exit_code == 0
		text = read(cli_outfile, String)
		@test Meta.parseall(text) isa Expr
		@test occursin("using Sunny", text)
		@test occursin("SpinWaveTheory", text)
	finally
		isfile(cli_outfile) && rm(cli_outfile)
	end
end

println("All MagestyCLI tests passed.")
