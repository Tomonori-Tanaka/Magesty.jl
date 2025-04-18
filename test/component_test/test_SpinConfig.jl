using .SpinConfigs
using LinearAlgebra

@testset "SpinConfig" begin
	# Test reading EMBSET file
	spinconfigs::Vector{SpinConfig} = read_embset("./examples/fept_tetragonal_2x2x2/EMBSET.dat", 16)

	# Test first configuration
	@test length(spinconfigs) >= 1  # At least 1 configurations in the file
	@test isapprox(spinconfigs[1].energy, -121.37424)
	@test isapprox(spinconfigs[2].energy, -121.37488)

	# Test first atom in first configuration
	@test isapprox(spinconfigs[1].magmom_size[1], norm([0.1870000, 0.0000000, 2.1670000]), atol=1e-4)
	@test isapprox(spinconfigs[1].spin_directions[:, 1], [0.1870000, 0.0000000, 2.1670000] / norm([0.1870000, 0.0000000, 2.1670000]), atol=1e-4)
	@test isapprox(spinconfigs[1].local_magfield[:, 1], [-1.17370e-02, 4.72560e-07, 1.01990e-03], atol=1e-4)

	# Test second atom in first configuration
	@test isapprox(spinconfigs[1].magmom_size[2], norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
	@test isapprox(spinconfigs[1].spin_directions[:, 2], [0.0000000, 0.0000000, 2.1760000] / norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
	@test isapprox(spinconfigs[1].local_magfield[:, 2], [-7.03430e-04, -5.49880e-06, -2.47420e-08], atol=1e-4)

	# Test ninth atom in first configuration (Pt atom)
	@test isapprox(spinconfigs[1].magmom_size[9], norm([0.0000000, 0.0000000, 0.2380000]), atol=1e-4)
	@test isapprox(spinconfigs[1].spin_directions[:, 9], [0.0000000, 0.0000000, 0.2380000] / norm([0.0000000, 0.0000000, 0.2380000]), atol=1e-4)
	@test isapprox(spinconfigs[1].local_magfield[:, 9], [9.35110e-03, 8.61050e-05, -4.37250e-06], atol=1e-4)

	# Test first atom in second configuration
	@test isapprox(spinconfigs[2].magmom_size[1], norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
	@test isapprox(spinconfigs[2].spin_directions[:, 1], [0.0000000, 0.0000000, 2.1760000] / norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
	@test isapprox(spinconfigs[2].local_magfield[:, 1], [1.01860e-03, -6.52810e-06, -5.18810e-08], atol=1e-4)

	# Test local_magfield_vertical calculation
	for config in spinconfigs
		for i in 1:size(config.local_magfield, 2)
			# Check that local_magfield_vertical is perpendicular to spin_directions
			@test isapprox(dot(config.local_magfield_vertical[:, i], config.spin_directions[:, i]), 0.0, atol=1e-10)
			
			# Check that local_magfield_vertical is correctly calculated
			proj = dot(config.local_magfield[:, i], config.spin_directions[:, i]) * config.spin_directions[:, i]
			expected_vertical = config.local_magfield[:, i] - proj
			@test isapprox(config.local_magfield_vertical[:, i], expected_vertical, atol=1e-10)
		end
	end

	# Test DataSet construction
	dataset = DataSet(spinconfigs, 0.8)
	@test dataset.training_data_num + dataset.validation_data_num == length(spinconfigs)
	@test length(dataset.training_data_indices) == dataset.training_data_num
	@test length(dataset.validation_data_indices) == dataset.validation_data_num
	@test all(i -> 1 ≤ i ≤ length(spinconfigs), dataset.training_data_indices)
	@test all(i -> 1 ≤ i ≤ length(spinconfigs), dataset.validation_data_indices)
	@test isempty(intersect(dataset.training_data_indices, dataset.validation_data_indices))

	# Test DataSet with training_ratio = 1.0
	dataset_all = DataSet(spinconfigs, 1.0)
	@test dataset_all.training_data_num == length(spinconfigs)
	@test isempty(dataset_all.validation_data_indices)
	@test dataset_all.training_data_indices == collect(1:length(spinconfigs))

	# Test error handling
	@test_throws ArgumentError DataSet(spinconfigs, 0.0)  # training_ratio <= 0
	@test_throws ArgumentError DataSet(spinconfigs, 1.1)  # training_ratio > 1
	@test_throws ArgumentError DataSet(SpinConfig[], 0.8)  # empty spinconfigs

	# Test parse_embset with specific indices
	use_indices = [1, 2]
	dataset_specific = parse_embset("./examples/fept_tetragonal_2x2x2/EMBSET.dat", 16, use_indices)
	@test length(dataset_specific.spinconfigs) == length(use_indices)
	@test all(i -> dataset_specific.spinconfigs[i].energy == spinconfigs[use_indices[i]].energy, 1:length(use_indices))

	# Test error handling for parse_embset
	@test_throws BoundsError parse_embset("./examples/fept_tetragonal_2x2x2/EMBSET.dat", 16, [0])  # index < 1
	@test_throws BoundsError parse_embset("./examples/fept_tetragonal_2x2x2/EMBSET.dat", 16, [length(spinconfigs)+1])  # index > length(spinconfigs)
end
