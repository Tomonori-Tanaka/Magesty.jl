"""
Regression tests for the `outcar_to_embset` API. The output must be
byte-identical to the golden EMBSET fixture, which was generated from the
OUTCAR fixture. A round trip through `read_embset` confirms the EMBSET text
parses back into spin configurations. The `magesty vasp embset` CLI
subcommand is covered separately by the `MagestyCLI` package tests.
"""

using Magesty
using Test

@testset "outcar_to_embset" begin
	dir         = joinpath(@__DIR__, "fixtures", "outcar")
	outcar_path = joinpath(dir, "OUTCAR")
	golden_text = read(joinpath(dir, "expected.embset"), String)

	# API function: the returned text is byte-identical to the golden file.
	@test outcar_to_embset([outcar_path]) == golden_text

	# The `output` keyword writes the same bytes to disk.
	api_outfile = tempname()
	try
		outcar_to_embset([outcar_path]; output = api_outfile)
		@test read(api_outfile, String) == golden_text
	finally
		isfile(api_outfile) && rm(api_outfile)
	end

	# Multiple OUTCARs produce one numbered block each; the EMBSET text
	# parses back through `read_embset`.
	multi_outfile = tempname()
	try
		outcar_to_embset([outcar_path, outcar_path]; output = multi_outfile)
		configs = read_embset(multi_outfile)
		@test length(configs) == 2
		@test configs[1].energy ≈ -10.5
		@test configs[2].energy ≈ -10.5
	finally
		isfile(multi_outfile) && rm(multi_outfile)
	end

	# An OUTCAR with no magnetic-moment block raises a clear error.
	no_magmom = tempname()
	try
		write(no_magmom, "no magnetic moment data here\n      1 F= -1.0 E0= -1.0\n")
		@test_throws ErrorException outcar_to_embset([no_magmom])
	finally
		isfile(no_magmom) && rm(no_magmom)
	end
end

println("All outcar_to_embset tests passed.")
