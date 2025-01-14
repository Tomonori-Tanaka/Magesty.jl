using TOML
using LinearAlgebra
@testset "b2feco1x1x2" begin
	input = """
	[general]
	name = "b2feco1x1x2"
	nat = 4
	kd = [ "Fe", "Co" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
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
		[ 0.0, 0.0, 2.0 ],
	]
	kd_list = [1, 1, 2, 2]
	position = [
		[ 0.00, 0.00, 0.00 ],
		[ 0.00, 0.00, 0.50 ],
		[ 0.50, 0.50, 0.25 ],
		[ 0.50, 0.50, 0.75 ], 
	]
	"""

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
	@test sclus.symmetry.international_symbol == "Pm-3m"
	@test sclus.symmetry.symnum_translation == [1, 17]
	@test sclus.symmetry.map_sym[1, 1] == 1
	@test sclus.symmetry.map_sym[1, 2] == 1
	@test sclus.symmetry.map_sym[3, 2] == 4
	@test sclus.symmetry.map_sym[3, 1] == 3
	@test sclus.symmetry.map_sym[1, 17] == 2
	@test sclus.symmetry.map_sym[2, 17] == 1
	@test sclus.symmetry.map_sym[3, 17] == 4
	@test sclus.symmetry.map_sym[4, 17] == 3


	# Magesty.BasisSets.__write_martix(sclus.basisset.projection_matrix)
	# idx = 17
	# @show sclus.symmetry.symdata[idx]
	# Magesty.Symmetries.__write_symdata(
	# 	sclus.symmetry.symdata,
	# 	"/Users/tomorin/Desktop",
	# 	"symdata.txt",
	# )
	# Magesty.BasisSets.__write_martix(sclus.basisset.each_projection_matrix[idx])
	# for (i, basis) in enumerate(sclus.basisset.basislist)
	# 	println(i, "\t", basis)
	# end

end



