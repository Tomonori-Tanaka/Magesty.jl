using TOML

@testset "l10feni1x1x2" begin
	input = """
			[general]
			name = "l10feni"
			nat = 4
			kd = [ "Fe", "Ni" ]
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
				Ni = [ 0, 1 ]
				[interaction.cutoff] # unit is bohr
				Fe-Fe = [ 0, -1 ] # first element is just dummy to align wigh lmax array
				Fe-Ni = [ 0, -1 ]
				Ni-Ni = [ 0, -1 ]
				# negative cutoff means all of the possible interaction will be considered.

			[regression]
			weight = 0.5        # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
			datafile = "EMBSET.txt"

			[structure]
			scale = 2.84
			lattice = [
				[ 1.0, 0.0, 0.0 ],
				[ 0.0, 1.0, 0.0 ],
				[ 0.0, 0.0, 3 ],
			]
			kd_list = [1, 1, 2, 2]
			position =[
				[0.00, 0.00, 0.00], 
				[0.00, 0.00, 0.50],
				[0.50, 0.50, 0.25],
				[0.50, 0.50, 0.75],
			]
			"""

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
	# for i in 1:sclus.symmetry.nsym
	# 	println(i)
    #     println(sclus.symmetry.symdata[i])
	# 	# display(sclus.symmetry.symdata[i].rotation_frac)
	# 	# println(sclus.symmetry.symdata[i].translation_frac)
	# end
end
