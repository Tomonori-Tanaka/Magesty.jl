using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/fept_tetragonal_2x2x2/input.toml", "r"),
)
input["regression"]["weight"] = 0.0
system = System(input, false)
sclus_energy = SpinCluster(system, false)
input["regression"]["weight"] = 1.0
sclus_torque = SpinCluster(system, false)

println(@sprintf("elapsed_time (structure): %10.6f", sclus_torque.structure.elapsed_time))
println(@sprintf("elapsed_time (symmetry):  %10.6f", sclus_torque.symmetry.elapsed_time))
println(@sprintf("elapsed_time (cluster):   %10.6f", sclus_torque.cluster.elapsed_time))
println(@sprintf("elapsed_time (basisset):  %10.6f", sclus_torque.basisset.elapsed_time))
println(@sprintf("elapsed_time (optimize):  %10.6f", sclus_torque.optimize.elapsed_time))
display(sclus_torque.optimize.SCE[5])
