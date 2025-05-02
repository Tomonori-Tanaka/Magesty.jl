function set_system(config)
	structure = Structure(
		config.lattice_vectors,
		config.is_periodic,
		config.kd_name,
		config.kd_int_list,
		config.x_fractional)
end

function set_symmetry(config, structure)
	symmetry = Symmetry(
		structure,
		config.tolerance_sym,
	)
end

function set_cluster(config, structure, symmetry)
	cluster = Cluster(structure, symmetry, config.nbody, config.cutoff_radii)
end

function set_basisset(config, structure, symmetry, cluster)
	basisset = BasisSet(structure, symmetry, cluster, config.lmax, config.nbody)
end

function set_optimize(config, structure, symmetry, basisset)
	optimize = SCEOptimizer(
		structure,
		symmetry,
		basisset,
		config.weight,
		config.datafile,
	)
end
