function set_system(parser)
	structure = Structure(
		parser.lattice_vectors,
		parser.is_periodic,
		parser.kd_name,
		parser.kd_int_list,
		parser.x_fractional)
end

function set_symmetry(parser, structure)
	symmetry = Symmetry(
		structure,
		parser.tolerance_sym,
	)
end

function set_cluster(parser, structure, symmetry)
	cluster = Cluster(structure, symmetry, parser.nbody, parser.cutoff_radii)
end

function set_basisset(parser, structure, symmetry, cluster)
	basisset = BasisSet(structure, symmetry, cluster, parser.lmax, parser.nbody)
end

function set_optimize(parser, structure, symmetry, basisset)
	optimize = SCEOptimizer(
		structure,
		symmetry,
		basisset,
		parser.j_zero_thr,
		parser.weight,
		parser.datafile,
	)
end
