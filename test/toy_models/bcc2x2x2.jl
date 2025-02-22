using Random
using LinearAlgebra
using TOML


@testset "toy_models/bcc2x2x2" begin
	input = """
	[general]
	name = "bcc2x2x2"
	mode = "check"
	nat = 16
	kd = [ "X" ]
	periodicity = [ true, true, true ]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-5

	[interaction]
		nbody = 2
		[interaction.lmax]
		X = [ 0, 1 ] # the number of elements shoud be the same with "nbody" value.
		[interaction.cutoff] # unit is bohr
		X-X = [ 0, 1.5 ] # first element is just dummy to align wigh lmax array
		# negative cutoff means all of the possible interaction will be considered.

	[regression]
	weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
	datafile = "EMBSET.txt"

	[structure]
	scale = 1
	lattice = [
		[ 2.0, 0.0, 0.0 ],
		[ 0.0, 2.0, 0.0 ],
		[ 0.0, 0.0, 2.0 ],
	]
	kd_list = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	position = [
		[ 0.00, 0.00, 0.00 ],
		[ 0.00, 0.00, 0.50 ],
		[ 0.00, 0.50, 0.00 ],
		[ 0.00, 0.50, 0.50 ],
		[ 0.50, 0.00, 0.00 ],
		[ 0.50, 0.00, 0.50 ],
		[ 0.50, 0.50, 0.00 ],
		[ 0.50, 0.50, 0.50 ],
		[ 0.25, 0.25, 0.25 ],
		[ 0.25, 0.25, 0.75 ],
		[ 0.25, 0.75, 0.25 ],
		[ 0.25, 0.75, 0.75 ],
		[ 0.75, 0.25, 0.25 ],
		[ 0.75, 0.25, 0.75 ],
		[ 0.75, 0.75, 0.25 ],
		[ 0.75, 0.75, 0.75 ],
	]
	"""
	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)

	# the number of atoms
	num_atoms = sclus.system.supercell.num_atoms
	# the number of independent basis functions
	num_basis = length(sclus.basisset.salc_list)
	# assign SCE values using random numbers
	min_val = -0.01 # -10 meV
	max_val = 0.01 # 10 meV
	sce_list = [min_val + rand() * (max_val - min_val) for _ in 1:num_basis+1]# including the bias term

	# generate configuration energy
	spinconfig_list = Vector{Magesty.SpinConfigs.SpinConfig}()
	num_config = 400
	design_matrix = zeros(num_config, num_basis + 1)
	design_matrix[:, 1] .= 1.0
	for i in 1:num_config
		# generate random spin configuration
		spin_config = rand(3, num_atoms)
		# normalize the spin moment
		for atom in 1:num_atoms
			norm_tmp = norm(spin_config[:, atom])
			spin_config[:, atom] /= norm_tmp
		end
		# calculate the design matrix
		for j in 1:num_basis
			design_matrix[i, j+1] = Magesty.Optimize.calc_design_matrix_element(
				sclus.basisset.salc_list[j],
				spin_config,
				sclus.symmetry,
			)
		end
		energy = dot(design_matrix[i, :], sce_list)
		push!(
			spinconfig_list,
			Magesty.SpinConfigs.SpinConfig(
				energy,
				ones(num_atoms),
				spin_config,
				zeros(3, num_atoms),
			),
		)
	end
	energy_list = design_matrix * sce_list

	# test the result
	parsed["general"]["mode"] = "optimize"
	optimize = Magesty.Optimize.SCEOptimizer(
		sclus.system,
		sclus.symmetry,
		sclus.basisset,
		sclus.config.j_zero_thr,
		sclus.config.weight,
		spinconfig_list,
	)
	sclus_new = SpinCluster(
		sclus.config,
		sclus.system,
		sclus.symmetry,
		sclus.cluster,
		sclus.basisset,
		optimize,
	)
	display(sclus_new.basisset.salc_list)
	# display(sclus_new.basisset.classified_basisdict)
	@test isapprox(sclus_new.optimize.energy_list, energy_list, atol=1e-6)
	@test isapprox(sclus_new.optimize.SCE, sce_list, atol=1e-6)
end
