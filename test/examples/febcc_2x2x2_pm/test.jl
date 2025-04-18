using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_pm/input.toml", "r"),
)
input["regression"]["weight"] = 0.0
system = System(input)
sclus_torque = SpinCluster(system)
input["regression"]["weight"] = 1.0
sclus_energy = SpinCluster(system, input, false)

display(sclus_torque.optimize.SCE)
display(sclus_energy.optimize.SCE)
