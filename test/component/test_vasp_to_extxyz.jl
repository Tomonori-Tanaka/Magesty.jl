"""
Regression tests for the `vasp_to_extxyz` API. The output must be
byte-identical to the golden extxyz fixtures, which were generated from
the VASP fixture inputs. The `magesty vasp extxyz` CLI subcommand is
covered separately by the `MagestyCLI` package tests.
"""

using Magesty
using Test

@testset "vasp_to_extxyz: $(system)" for system in ("FeRh", "IrMn3")
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
end

println("All vasp_to_extxyz tests passed.")
