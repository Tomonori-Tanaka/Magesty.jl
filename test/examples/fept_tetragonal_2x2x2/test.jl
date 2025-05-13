using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/fept_tetragonal_2x2x2/input.toml", "r"),
)
system = System(input, false)
sclus = SpinCluster(system, input, false)
Magesty.write_sce2xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))

# weight_list = collect(0.1:0.1:1.0)
# for weight in weight_list
# 	input["regression"]["weight"] = weight
# 	sclus_restart = SpinCluster(sclus, input, false)
# 	println("weight: ", weight)
# 	println("1NN: ", sclus_restart.optimize.SCE[2])	
# 	println("relative error energy:            ", sclus_restart.optimize.relative_error_energy)
# 	println("relative error magfield_vertical: ", sclus_restart.optimize.relative_error_magfield_vertical)
# end
# println(@sprintf("elapsed_time (structure): %10.6f", sclus_torque.structure.elapsed_time))
# println(@sprintf("elapsed_time (symmetry):  %10.6f", sclus_torque.symmetry.elapsed_time))
# println(@sprintf("elapsed_time (cluster):   %10.6f", sclus_torque.cluster.elapsed_time))
# println(@sprintf("elapsed_time (basisset):  %10.6f", sclus_torque.basisset.elapsed_time))
# println(@sprintf("elapsed_time (optimize):  %10.6f", sclus_torque.optimize.elapsed_time))
