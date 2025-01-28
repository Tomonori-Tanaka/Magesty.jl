using TOML

@testset "b20fege1x1x1" begin
	input = """
	[general]
	name = "b20fege"
	nat = 8
	kd = [ "Fe", "Ge" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-5

	[interaction]
		nbody = 2
		[interaction.lmax]
		Fe = [ 0, 1 ] # the number of elements shoud be the same with "nbody" value.
		Ge = [ 0, 0 ]
		[interaction.cutoff] # unit is bohr
		Fe-Fe = [ 0, -1 ] # first element is just dummy to align wigh lmax array
		Fe-Ge = [ 0, -1 ]
		Ge-Ge = [ 0, 0 ]
		# negative cutoff means all of the possible interaction will be considered.

	[regression]
	weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
	datafile = "EMBSET.txt"

	[structure]
	scale = 4.6647157669
	lattice = [
		[ 1.0, 0.0, 0.0 ],
		[ 0.0, 1.0, 0.0 ],
		[ 0.0, 0.0, 1.0 ],
	]
	kd_list = [
	1, 1, 1, 1, 
	2, 2, 2, 2,
	]
	position = [
	[0.135711998,        0.135711998,        0.135711998],
	[0.364288002,        0.864287972,        0.635712028],
	[0.635712028,        0.364288002,        0.864287972],
	[0.864287972,        0.635712028,        0.364288002],
	[0.842083991,        0.842083991,        0.842083991],
	[0.657916009,        0.157915995,        0.342083991],
	[0.342083991,        0.657916009,        0.157915995],
	[0.157915995,        0.342083991,        0.657916009],
]
	"""
	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
end
