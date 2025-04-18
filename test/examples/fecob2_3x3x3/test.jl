using TOML

input = TOML.parse(open("/Users/tomorin/Packages/Magesty/test/examples/fecob2_3x3x3/input.toml", "r"))
system = System(input)
sclus = SpinCluster(system)
display(sclus.optimize.SCE)
input["regression"]["weight"] = 1.0
system_torque = System(input)
sclus_torque = SpinCluster(system_torque)
display(sclus_torque.optimize.SCE)