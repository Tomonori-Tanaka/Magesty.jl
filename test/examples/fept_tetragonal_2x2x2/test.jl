using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/fept_tetragonal_2x2x2/input.toml", "r"),
)
system = System(input, verbosity = false)
sclus = SpinCluster(system, input, verbosity = false)
Magesty.write_sce2xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))
