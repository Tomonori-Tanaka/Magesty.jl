using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_fm/input.toml", "r"),
)
system = System(input)
sclus = SpinCluster(system)
print_info(sclus)
display(sclus.optimize.SCE)