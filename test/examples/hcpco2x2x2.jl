using TOML

@testset "hcpco2x2x2" begin
	input = """
	[general]
	name = "hcpco"
	nat = 16
	kd = [ "Co" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-4

	[interaction]
		nbody = 2
		[interaction.lmax]
		Co = [ 0, 1 ]
		[interaction.cutoff] # unit is bohr
		Co-Co = [ 0, 2.48 ]
		# negative cutoff means all of the possible interaction will be considered.

	[regression]
	weight = 0.5            # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
	datafile = "EMBSET.txt"

	[structure]
	scale = 1.0
	lattice = [
	[4.9774985313,        0.0000000000,        0.0000000000],
	[2.4887495527,        4.3106400097,        0.0000000000],
	[0.0000000000,        0.0000000000,        8.0566759109],
	]
	kd_list = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	position = [
	[0.000000000,        0.000000000,        0.000000000],
	[0.000000000,        0.000000000,        0.500000000],
	[0.000000000,        0.500000000,        0.000000000],
	[0.000000000,        0.500000000,        0.500000000],
	[0.500000000,        0.000000000,        0.000000000],
	[0.500000000,        0.000000000,        0.500000000],
	[0.500000000,        0.500000000,        0.000000000],
	[0.500000000,        0.500000000,        0.500000000],
	[0.166666672,        0.166666672,        0.250000000],
	[0.166666672,        0.166666672,        0.750000000],
	[0.166666672,        0.666666687,        0.250000000],
	[0.166666672,        0.666666687,        0.750000000],
	[0.666666687,        0.166666672,        0.250000000],
	[0.666666687,        0.166666672,        0.750000000],
	[0.666666687,        0.666666687,        0.250000000],
	[0.666666687,        0.666666687,        0.750000000]
	]
	"""
	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)

	# for data in sclus.symmetry.symdata
	# 	@show data
	# end
	# for itrans in sclus.symmetry.symnum_translation
	# 	@show sclus.symmetry.symdata[itrans].translation_frac
	# end
end
