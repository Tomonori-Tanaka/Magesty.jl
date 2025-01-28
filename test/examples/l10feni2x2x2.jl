using TOML

@testset "l10feni2x2x2" begin
	input = """
	[general]
	name = "l10feni2x2x2"
	nat = 16
	kd = [ "Fe", "Ni" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-5

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
	scale = 1
	lattice = [
		[ 2.0, 0.0, 0.0 ],
		[ 0.0, 2.0, 0.0 ],
		[ 0.0, 0.0, 2.2 ],
	]
	kd_list = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2]
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
	# for i in 1:sclus.symmetry.nsym
	# 	if sclus.symmetry.symdata[i].is_translation_included
	# 		continue
	# 	end
	# 	println(i)
	# 	println(sclus.symmetry.symdata[i])
	# end

	# idx = 48
	# Magesty.BasisSets.__write_martix(sclus.basisset.each_projection_dict[1][idx])
	# for idx in 1:sclus.symmetry.nsym
	# 	display(sclus.basisset.each_projection_dict[1][idx])
	# end
	# for (i, basis) in enumerate(sclus.basisset.basislist)
	# 	println(i, "\t", basis)
	# println(sclus.basisset.basislist)
	# @test length(sclus.basisset.basislist) ==
	# 	  8 * (3 * 3) + (8 - 1) * (3 * 3) + 2 * (3 * (3 * 3)) + 2 * (4 * (3 * 3)) + 2 # 72 + 63 + 54 +72 + 2 
end
