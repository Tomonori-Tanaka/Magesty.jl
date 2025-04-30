using TOML
using Printf

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_pm/input.toml", "r"),
)
input["regression"]["weight"] = 0.0
system = System(input, false)
sclus_torque = SpinCluster(system, false)
input["regression"]["weight"] = 1.0
sclus_energy = SpinCluster(system, input, false)

println("1NN (torque): ", sclus_torque.optimize.SCE[2])
println("1NN (energy): ", sclus_energy.optimize.SCE[2])
println(@sprintf("elapsed_time (structure): %10.6f", sclus_torque.structure.elapsed_time))
println(@sprintf("elapsed_time (symmetry):  %10.6f", sclus_torque.symmetry.elapsed_time))
println(@sprintf("elapsed_time (cluster):   %10.6f", sclus_torque.cluster.elapsed_time))
println(@sprintf("elapsed_time (basisset):  %10.6f", sclus_torque.basisset.elapsed_time))
println(@sprintf("elapsed_time (optimize):  %10.6f", sclus_torque.optimize.elapsed_time))
