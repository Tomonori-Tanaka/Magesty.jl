using TOML

@testset "hcpco1x1x1" begin
	input = """
	[general]
	name = "hcpco"
	nat = 2
	kd = [ "Co" ]
	periodicity = [ true, true, true ]
	# periodicity = [false, false, false]
	j_zero_thr = 1e-10

	[symmetry]
	tolerance = 1e-8
	time-reversal = true
	crystal = true

	[interaction]
		nbody = 2
		[interaction.lmax]
		Co = [ 2, 1 ]
		[interaction.cutoff] # unit is bohr
		Co-Co = [ 0, -1 ]
		# negative cutoff means all of the possible interaction will be considered.

	[regression]
	weight = 0.5            # 1.0: only use the energy info. 0.0: only use torque info (0 <= weight_ratio <= 1)
	datafile = "EMBSET.txt"

	[structure]
	scale = 1.0
	lattice = [
		[ 1.24437467257, 2.15532015653, 0.0000000000 ],
		[ -1.24437467257, 2.15532015653, 0.0000000000 ],
		[ 0.0000000000, 0.0000000000, 4.02833795545 ],
	]
        kd_list = [1, 1]
	position = [
		[ 0.000000000, 0.000000000, 0.000000000 ],
		[ 0.333333333, 0.333333333, 0.500000000 ],
	]
	"""
	parsed = TOML.parse(input)
	sclus = SpinCluster(parsed)
end
