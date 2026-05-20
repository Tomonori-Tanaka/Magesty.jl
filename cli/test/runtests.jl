"""
Regression test for the `magesty vasp extxyz` CLI subcommand. The command
output must be byte-identical to the golden extxyz fixtures, which live in
the core package test tree (`test/component/fixtures/`) and are reused
here by resolving a path relative to the repository root.
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

println("All MagestyCLI tests passed.")
