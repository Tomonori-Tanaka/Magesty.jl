using TOML

input = TOML.parse(open("/Users/tomorin/Packages/Magesty/test/examples/feptL10_2x2x2_saxis001/input.toml", "r"))
sclus = SpinCluster(input)