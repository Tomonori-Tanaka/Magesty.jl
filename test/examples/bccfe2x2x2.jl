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
   model = 1 # 1: isotropic Heisenberg model
	   # nbody = 2
	   [interaction.lmax]
	   Fe = [ 0, 1 ] # the number of elements shoud be the same with "nbody" value.
	   [interaction.cutoff] # unit is bohr
	   Fe-Fe = [ 0.0, 2.8 ] # first element is just dummy to align wigh lmax array

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
   position = [
	   { index = 1, kd = 1, coords = [ 0.00, 0.00, 0.00 ] },
	   { index = 2, kd = 1, coords = [ 0.25, 0.25, 0.25 ] },
	   { index = 3, kd = 1, coords = [ 0.00, 0.00, 0.50 ] },
	   { index = 4, kd = 1, coords = [ 0.50, 0.00, 0.00 ] },
	   { index = 5, kd = 1, coords = [ 0.00, 0.50, 0.00 ] },
	   { index = 6, kd = 1, coords = [ 0.25, 0.25, 0.75 ] },
	   { index = 7, kd = 1, coords = [ 0.75, 0.25, 0.25 ] },
	   { index = 8, kd = 1, coords = [ 0.25, 0.75, 0.25 ] },
	   { index = 9, kd = 1, coords = [ 0.00, 0.50, 0.50 ] },
	   { index = 10, kd = 1, coords = [ 0.50, 0.00, 0.50 ] },
	   { index = 11, kd = 1, coords = [ 0.50, 0.50, 0.00 ] },
	   { index = 12, kd = 1, coords = [ 0.25, 0.75, 0.75 ] },
	   { index = 13, kd = 1, coords = [ 0.75, 0.25, 0.75 ] },
	   { index = 14, kd = 1, coords = [ 0.75, 0.75, 0.25 ] },
	   { index = 15, kd = 1, coords = [ 0.50, 0.50, 0.50 ] },
	   { index = 16, kd = 1, coords = [ 0.75, 0.75, 0.75 ] },
   ]
   """

	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
	println(sclus.basisset.basislist)


end