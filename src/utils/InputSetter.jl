function set_system(parser)
	system = System(
		parser.lattice_vectors,
		parser.is_periodic,
		parser.kd_name,
		parser.kd_int_list,
		parser.x_fractional)
end

function set_symmetry(parser, system)
	symmetry = Symmetry(
		system,
		parser.tolerance_sym,
	)
end

function set_cluster(parser, system, symmetry)
	cluster = Cluster(system, symmetry, parser.nbody, parser.cutoff_radii)
end

function set_basisset(parser, system, symmetry, cluster)
	basisset = BasisSet(system, symmetry, cluster, parser.lmax, parser.nbody)
end

function set_optimize(parser, system, symmetry, cluster)
	optimize = Optimize(
		system,
		symmetry,
		cluster,
		parser.time_reversal_sym,
		parser.crystal_sym,
		parser.lmax,
		parser.model,
	)
end
