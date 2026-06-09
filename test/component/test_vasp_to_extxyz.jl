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

# A non-converged run (scstep count >= NELM) is still written, tagged
# `converged=F`, and warns. We synthesize one by lowering NELM below the
# fixture's scstep count, so the count rule reports non-convergence.
@testset "vasp_to_extxyz: non-converged run warns and tags F" begin
    src = read(joinpath(@__DIR__, "fixtures", "FeRh", "vasprun.xml"), String)
    # FeRh runs 35 electronic steps; force the cap to 30 so 35 >= 30 → F.
    edited = replace(src, "name=\"NELM\">   100" => "name=\"NELM\">    30")
    @test edited != src   # guard: the substitution actually fired

    tmp = tempname() * ".xml"
    try
        write(tmp, edited)
        text = @test_logs (:warn, r"did not converge") vasp_to_extxyz(tmp)
        header = split(text, "\n")[2]
        @test occursin("converged=F", header)
    finally
        isfile(tmp) && rm(tmp)
    end
end

println("All vasp_to_extxyz tests passed.")
