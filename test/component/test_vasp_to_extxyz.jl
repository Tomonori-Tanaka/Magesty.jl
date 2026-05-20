"""
Regression tests for `vasp_to_extxyz` and the `magesty vasp extxyz`
subcommand. The output must be byte-identical to the golden extxyz
fixtures, which were generated from the VASP fixture inputs.
"""

using Magesty
using Test

@testset "vasp_to_extxyz / vasp extxyz: $(system)" for system in ("FeRh", "IrMn3")
	dir          = joinpath(@__DIR__, "fixtures", system)
	vasprun_path = joinpath(dir, "vasprun.xml")
	oszicar_path = joinpath(dir, "OSZICAR")
	golden_text  = read(joinpath(dir, "$(system).extxyz"), String)

	# API function: the returned text is byte-identical to the golden file.
	@test vasp_to_extxyz(vasprun_path; oszicar = oszicar_path) == golden_text

	# The `output` keyword writes the same bytes to disk.
	api_outfile = tempname() * ".extxyz"
	try
		vasp_to_extxyz(vasprun_path; oszicar = oszicar_path, output = api_outfile)
		@test read(api_outfile, String) == golden_text
	finally
		isfile(api_outfile) && rm(api_outfile)
	end

	# CLI subcommand `magesty vasp extxyz`: exits 0 and writes the same bytes.
	cli_outfile = tempname() * ".extxyz"
	try
		exit_code = Magesty.command_main(
			["vasp", "extxyz", vasprun_path, "--oszicar", oszicar_path, "--output", cli_outfile],
		)
		@test exit_code == 0
		@test read(cli_outfile, String) == golden_text
	finally
		isfile(cli_outfile) && rm(cli_outfile)
	end
end

println("All vasp_to_extxyz tests passed.")
