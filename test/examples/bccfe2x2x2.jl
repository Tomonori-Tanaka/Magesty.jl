using TOML
@testset "bccfe2x2x2" begin
	input = """
	[general]
	name = "bccfe"
	nat = 16
	kd = [ "Fe" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-8

	[interaction]
	nbody = 2
		[interaction.lmax]
		Fe = [ 0, 1 ] # the number of elements shoud be the same with "nbody" value.
		[interaction.cutoff] # unit is bohr
	  	Fe-Fe = [ 0.0, 2.80 ] # first element is just dummy to align wigh lmax array

	[regression]
	weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
	datafile = "EMBSET.txt"

	[structure]
	scale = 2.83
	lattice = [
		[ 2.0, 0.0, 0.0 ],
		[ 0.0, 2.0, 0.0 ],
		[ 0.0, 0.0, 2.0 ],
	]
	kd_list = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	position = [
		[ 0.00, 0.00, 0.00 ],
		[ 0.25, 0.25, 0.25 ],
		[ 0.00, 0.00, 0.50 ],
	 	[ 0.50, 0.00, 0.00 ],
		[ 0.00, 0.50, 0.00 ],
		[ 0.25, 0.25, 0.75 ],
		[ 0.75, 0.25, 0.25 ],
		[ 0.25, 0.75, 0.25 ],
		[ 0.00, 0.50, 0.50 ],
		[ 0.50, 0.00, 0.50 ],
		[ 0.50, 0.50, 0.00 ],
		[ 0.25, 0.75, 0.75 ],
		[ 0.75, 0.25, 0.75 ],
		[ 0.75, 0.75, 0.25 ],
		[ 0.50, 0.50, 0.50 ],
		[ 0.75, 0.75, 0.75 ],
	]
   """

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
	# println(sclus.basisset.basislist)


end