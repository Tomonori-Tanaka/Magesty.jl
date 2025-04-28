using TOML

input = TOML.parse(open("/Users/tomorin/Packages/Magesty/test/examples/fecob2_3x3x3/input.toml", "r"))
# system = System(input)
# sclus = SpinCluster(system, false)
# display(sclus.optimize.SCE)
input["regression"]["weight"] = 1.0
system_torque = System(input, false)
sclus_torque = SpinCluster(system_torque, false)
println("1NN (torque): ", sclus_torque.optimize.SCE[2])
println(@sprintf("elapsed_time (structure): %10.6f", sclus_torque.structure.elapsed_time))
println(@sprintf("elapsed_time (symmetry):  %10.6f", sclus_torque.symmetry.elapsed_time))
println(@sprintf("elapsed_time (cluster):   %10.6f", sclus_torque.cluster.elapsed_time))
println(@sprintf("elapsed_time (basisset):  %10.6f", sclus_torque.basisset.elapsed_time))
println(@sprintf("elapsed_time (optimize):  %10.6f", sclus_torque.optimize.elapsed_time))
