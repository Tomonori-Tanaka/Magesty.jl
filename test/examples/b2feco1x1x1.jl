using TOML

@testset "b2feco1x1x1" begin
	input = """
			[general]
			name = "b2feco"
			nat = 2
			kd = [ "Fe", "Co" ]
			periodicity = [ true, true, true ]
			# periodicity = [false, false, false]
			# periodicity = [true, false, false]
			j_zero_thr = 1e-10

			[symmetry]
			tolerance = 1e-8

			[interaction]
			nbody = 2
				[interaction.lmax]
				Fe = [ 0, 1 ] # the number of elements shoud be the same with "nbody" value.
				Co = [ 0, 1 ]
				[interaction.cutoff] # unit is bohr
				Fe-Fe = [ 0, -1 ] # first element is just dummy to align wigh lmax array
				Fe-Co = [ 0, -1 ]
				Co-Co = [ 0, -1 ]
				# negative cutoff means all of the possible interaction will be considered.

			[regression]
			weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
			datafile = "EMBSET.txt"

			[structure]
			scale = 2.84
			lattice = [
				[ 1.0, 0.0, 0.0 ],
				[ 0.0, 1.0, 0.0 ],
				[ 0.0, 0.0, 1.0 ],
			]
			kd_list = [1, 2]
			position =[
				[ 0.00, 0.00, 0.00 ], 
				[ 0.50, 0.50, 0.50 ],
			]
			"""

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
	println(sclus.basisset.basislist)

end
