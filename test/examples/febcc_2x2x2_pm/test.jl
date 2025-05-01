using TOML
using Printf

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_pm/input.toml", "r"),
)
system = System(input, false)
weight_list = collect(0.0:0.2:1.0)
for weight in weight_list
	println("weight: ", weight)
	input["regression"]["weight"] = weight
	sclus = SpinCluster(system, input, false)
	println("1NN: ", sclus.optimize.SCE[2])	
	println("relative error energy:            ", sclus.optimize.relative_error_energy)
	println("relative error magfield_vertical: ", sclus.optimize.relative_error_magfield_vertical)
end

# println(@sprintf("elapsed_time (structure): %10.6f", sclus_torque.structure.elapsed_time))
# println(@sprintf("elapsed_time (symmetry):  %10.6f", sclus_torque.symmetry.elapsed_time))
# println(@sprintf("elapsed_time (cluster):   %10.6f", sclus_torque.cluster.elapsed_time))
# println(@sprintf("elapsed_time (basisset):  %10.6f", sclus_torque.basisset.elapsed_time))
# println(@sprintf("elapsed_time (optimize):  %10.6f", sclus_torque.optimize.elapsed_time))
