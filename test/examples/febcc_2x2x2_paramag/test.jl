using TOML

input = TOML.parse(open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_paramag/input.toml", "r"))
sclus = SpinCluster(input)