using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_fm/input.toml", "r"),
)
input["regression"]["weight"] = 0.0
system_energy = System(input)
sclus_energy = SpinCluster(system_energy)
input["regression"]["weight"] = 1.0
system_torque = System(input)
sclus_torque = SpinCluster(system_torque)

display(sclus_energy.optimize.SCE)
display(sclus_torque.optimize.SCE)
