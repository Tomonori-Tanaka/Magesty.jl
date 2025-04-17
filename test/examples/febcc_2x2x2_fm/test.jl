using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_fm/input.toml", "r"),
)
sclus = SpinCluster(input)
print_info(sclus)
display(sclus.optimize.SCE)