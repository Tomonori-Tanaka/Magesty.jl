using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_fm/input.toml", "r"),
)
system_energy = System(input)
sclus = SpinCluster(system_energy)
print_info(sclus)
display(sclus.optimize.SCE)
input["regression"]["weight"] = 1.0
system_torque = System(input)
sclus_torque = SpinCluster(system_torque)
display(sclus_torque.optimize.SCE)
