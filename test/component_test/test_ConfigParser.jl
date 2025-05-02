using .ConfigParser
using Test

@testset "ConfigParser" begin
	@testset "Config4System" begin
		# Create a valid input dictionary
		input_dict = Dict{String, Any}(
			"general" => Dict{String, Any}(
				"name" => "test_system",
				"nat" => 2,
				"kd" => ["H", "O"],
				"periodicity" => [true, true, true],
			),
			"symmetry" => Dict{String, Any}(
				"tolerance" => 1e-3,
			),
			"interaction" => Dict{String, Any}(
				"nbody" => 2,
				"lmax" => Dict{String, Any}(
					"H" => [1, 2],
					"O" => [3, 4],
				),
				"cutoff" => Dict{String, Any}(
					"H-H" => [0.0, 2.0],
					"H-O" => [0.0, 3.0],
					"O-O" => [0.0, 4.0],
				),
			),
			"structure" => Dict{String, Any}(
				"lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
				"kd_list" => [1, 2],
				"position" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
			),
		)

		# Test valid input
		config = Config4System(input_dict)
		@test config.name == "test_system"
		@test config.num_atoms == 2
		@test config.kd_name == ["H", "O"]
		@test config.nbody == 2
		@test config.lmax == [1 2; 3 4]
		@test config.cutoff_radii[1, 1, 2] ≈ 2.0
		@test config.cutoff_radii[1, 2, 2] ≈ 3.0
		@test config.cutoff_radii[2, 2, 2] ≈ 4.0
		@test config.lattice_vectors == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
		@test config.kd_int_list == [1, 2]
		@test config.x_fractional[:, 1] ≈ [0.0, 0.0, 0.0]
		@test config.x_fractional[:, 2] ≈ [0.5, 0.5, 0.5]
		@test config.is_periodic == [true, true, true]
		@test config.tolerance_sym ≈ 1e-3

		# Test missing required section
		invalid_dict = copy(input_dict)
		delete!(invalid_dict, "general")
		@test_throws ArgumentError Config4System(invalid_dict)

		# Test invalid name
		invalid_dict = copy(input_dict)
		invalid_dict["general"]["name"] = ""
		@test_throws ArgumentError Config4System(invalid_dict)

		# Test invalid number of atoms
		invalid_dict = copy(input_dict)
		invalid_dict["general"]["nat"] = 0
		@test_throws ArgumentError Config4System(invalid_dict)
	end

	@testset "Config4Optimize" begin
		# Create a valid input dictionary
		input_dict = Dict{String, Any}(
			"regression" => Dict{String, Any}(
				"datafile" => "test.dat",
				"ndata" => 100,
				"weight" => 0.5,
			),
		)

		# Test valid input
		config = Config4Optimize(input_dict)
		@test config.datafile == "test.dat"
		@test config.ndata == 100
		@test config.weight ≈ 0.5

		# Test missing required section
		invalid_dict = copy(input_dict)
		delete!(invalid_dict, "regression")
		@test_throws ArgumentError Config4Optimize(invalid_dict)

		# Test empty datafile name
		invalid_dict = copy(input_dict)
		invalid_dict["regression"]["datafile"] = ""
		@test_throws ArgumentError Config4Optimize(invalid_dict)

		# Test invalid weight
		invalid_dict = copy(input_dict)
		invalid_dict["regression"]["weight"] = 1.5
		@test_throws ArgumentError Config4Optimize(invalid_dict)
	end

	@testset "parse_lmax" begin
		kd_name = ["H", "O"]
		nbody = 2
		lmax_dict = Dict{String, Any}(
			"H" => [1, 2],
			"O" => [3, 4],
		)

		lmax = ConfigParser.parse_lmax(lmax_dict, kd_name, nbody)
		@test size(lmax) == (2, 2)
		@test lmax[1, 1] == 1
		@test lmax[1, 2] == 2
		@test lmax[2, 1] == 3
		@test lmax[2, 2] == 4

		# Test error cases
		invalid_dict = Dict{String, Any}("H" => [1])
		@test_throws ArgumentError ConfigParser.parse_lmax(invalid_dict, kd_name, nbody)
	end

	@testset "parse_cutoff" begin
		kd_name = ["H", "O"]
		nbody = 2
		cutoff_dict = Dict{String, Any}(
			"H-H" => [0.0, 2.0],
			"H-O" => [0.0, 3.0],
			"O-O" => [0.0, 4.0],
		)

		cutoff = ConfigParser.parse_cutoff(cutoff_dict, kd_name, nbody)
		@test size(cutoff) == (2, 2, 2)
		@test cutoff[1, 1, 2] ≈ 2.0
		@test cutoff[1, 2, 2] ≈ 3.0
		@test cutoff[2, 2, 2] ≈ 4.0
		@test cutoff[2, 1, 2] ≈ 3.0  # Check symmetry

		# Test error cases
		invalid_dict = Dict{String, Any}("H-H" => [0.0])
		@test_throws ArgumentError ConfigParser.parse_cutoff(invalid_dict, kd_name, nbody)
	end

	@testset "parse_position" begin
		num_atoms = 2
		position_list = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

		positions = ConfigParser.parse_position(position_list, num_atoms)
		@test size(positions) == (3, 2)
		@test positions[:, 1] ≈ [0.0, 0.0, 0.0]
		@test positions[:, 2] ≈ [0.5, 0.5, 0.5]

		# Test error cases
		invalid_list = [[0.0, 0.0]]
		@test_throws ArgumentError ConfigParser.parse_position(invalid_list, num_atoms)
	end
end
