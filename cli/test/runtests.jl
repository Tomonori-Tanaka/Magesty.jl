"""
Regression tests for the `magesty vasp` CLI subcommands (`extxyz`, `toml`,
`embset`). The command output must be byte-identical to the golden
fixtures, which live in the core package test tree
(`test/component/fixtures/`) and are reused here by resolving a path
relative to the repository root.
"""

using MagestyCLI
using Magesty: IncarIO
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
	dir          = joinpath(FIXTURE_DIR, "oszicar")
	oszicar_path = joinpath(dir, "OSZICAR")
	golden_text  = read(joinpath(dir, "expected.embset"), String)

	# CLI subcommand `magesty vasp embset`: exits 0 and writes the same bytes.
	cli_outfile = tempname()
	try
		exit_code = MagestyCLI.command_main(
			["vasp", "embset", oszicar_path, "--output", cli_outfile],
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

@testset "magesty vasp mfa" begin
	incar_path = joinpath(FIXTURE_DIR, "incar", "INCAR")

	# CLI subcommand `magesty vasp mfa`: exits 0 and writes
	# num_points * num_samples INCAR files. Sampling is stochastic, so we
	# assert structural invariants (file count, naming, magnitude
	# preservation, MAGMOM == M_CONSTR) rather than byte-exact output.
	mktempdir() do dir
		outdir = joinpath(dir, "out")
		exit_code = MagestyCLI.command_main(
			["vasp", "mfa", incar_path, "tau",
				"--start", "0.1", "--stop", "0.3", "--num-points", "3",
				"--num-samples", "2", "--outdir", outdir],
		)
		@test exit_code == 0
		files = sort(readdir(outdir))
		@test files == ["sample-$(i).INCAR" for i = 1:6]
		for f in files
			d = IncarIO.parse_incar(joinpath(outdir, f))
			m = d[:MAGMOM]
			@test d[:M_CONSTR] == m
			for i = 1:4
				mag = sqrt(sum(abs2, m[3i-2:3i]))
				@test isapprox(mag, 3.0; atol = 1e-6)
			end
		end
	end
end

println("All MagestyCLI tests passed.")
