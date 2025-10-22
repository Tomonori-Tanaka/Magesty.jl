using .ConfigParser
using Test

@testset "ConfigParser" begin
	@testset "Config4System basic and interaction parsing" begin
		# Create a valid input dictionary matching latest schema
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
				"body1" => Dict{String, Any}(
					"lmax" => Dict{String, Any}("H" => 1, "O" => 3),
				),
				"body2" => Dict{String, Any}(
					"lsum" => 4,
					"cutoff" => Dict{String, Any}(
						"H-H" => 2.0,
						"H-O" => 3.0,
						"O-O" => 4.0,
					),
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
		@test config.body1_lmax == [1, 3]
		# bodyn_lsum is OffsetArray indexed from 2
		@test config.bodyn_lsum[2] == 4
		# bodyn_cutoff is symmetric and indexed as [n, i, j]
		@test config.bodyn_cutoff[2, 1, 1] ≈ 2.0
		@test config.bodyn_cutoff[2, 1, 2] ≈ 3.0
		@test config.bodyn_cutoff[2, 2, 1] ≈ 3.0
		@test config.bodyn_cutoff[2, 2, 2] ≈ 4.0
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

		# Empty name is currently allowed (no validation in constructor)
		invalid_dict = copy(input_dict)
		invalid_dict["general"]["name"] = ""
		config_empty = Config4System(invalid_dict)
		@test config_empty.name == ""

		# Test invalid number of atoms
		invalid_dict = copy(input_dict)
		invalid_dict["general"]["nat"] = 0
		@test_throws ArgumentError Config4System(invalid_dict)
	end

	@testset "Config4System nbody=1 edge case" begin
		input_dict = Dict{String, Any}(
			"general" => Dict{String, Any}(
				"name" => "nbody1",
				"nat" => 1,
				"kd" => ["Fe"],
			),
			"symmetry" => Dict{String, Any}(),
			"interaction" => Dict{String, Any}(
				"nbody" => 1,
				"body1" => Dict{String, Any}("lmax" => Dict{String, Any}("Fe" => 2)),
			),
			"structure" => Dict{String, Any}(
				"lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
				"kd_list" => [1],
				"position" => [[0.0, 0.0, 0.0]],
			),
		)
		config = Config4System(input_dict)
		@test config.nbody == 1
		@test config.body1_lmax == [2]
		# bodyn_* should be empty ranges for nbody=1
		@test length(config.bodyn_lsum) == 0
		@test size(parent(config.bodyn_cutoff)) == (0, 1, 1)
	end

	@testset "Config4Optimize" begin
		# Create a valid input dictionary
		input_dict = Dict{String, Any}(
			"regression" => Dict{String, Any}(
				"datafile" => "test.dat",
				"ndata" => 100,
				"weight" => 0.5,
				"alpha" => 0.1,
				"lambda" => 0.01,
			),
		)

		# Test valid input
		config = Config4Optimize(input_dict)
		@test config.datafile == "test.dat"
		@test config.ndata == 100
		@test config.weight ≈ 0.5
		@test config.alpha ≈ 0.1
		@test config.lambda ≈ 0.01

		# Test missing required section
		invalid_dict = copy(input_dict)
		delete!(invalid_dict, "regression")
		@test_throws ArgumentError Config4Optimize(invalid_dict)

		# Test empty datafile name
		invalid_dict = Dict{String, Any}("regression" => Dict{String, Any}("datafile" => ""))
		@test_throws ArgumentError Config4Optimize(invalid_dict)

		# Test invalid weight
		invalid_dict = Dict{String, Any}(
			"regression" => Dict{String, Any}(
				"datafile" => "test.dat",
				"weight" => 1.5,
			),
		)
		@test_throws ArgumentError Config4Optimize(invalid_dict)
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
