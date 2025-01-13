using TOML

@testset "fesi_atomicchain" begin
	input = """
			[general]
			name = "fesi_atomicchain"
			nat = 4
			kd = [ "Fe", "Si" ]
			# periodicity = [ true, true, true ]
			# periodicity = [false, false, false]
			periodicity = [true, false, false]
			j_zero_thr = 1e-10

			[symmetry]
			tolerance = 1e-8

			[interaction]
			nbody = 2
				[interaction.lmax]
				Fe = [ 2, 1 ] # the number of elements shoud be the same with "nbody" value.
				Si = [ 2, 0 ]
				[interaction.cutoff] # unit is bohr
				Fe-Fe = [ 0, -1 ] # first element is just dummy to align wigh lmax array
				Fe-Si = [ 0, -1 ]
				Si-Si = [ 0, -1 ]
				# negative cutoff means all of the possible interaction will be considered.

			[regression]
			weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
			datafile = "EMBSET.txt"

			[structure]
			scale = 1
			lattice = [
				[ 4.0, 0.0, 0.0 ],
				[ 0.0, 10.0, 0.0 ],
				[ 0.0, 0.0, 10.0 ],
			]
			kd_list = [1, 1, 2, 2]
			position =[
				[ 0.00, 0.00, 0.10 ], 
                [ 0.40, 0.00, 0.10 ],
				[ 0.00, 0.00, 0.00 ],
                [ 0.40, 0.00, 0.00 ],
			]
			"""

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)

end