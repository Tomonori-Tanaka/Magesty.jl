using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty.Optimize
using Magesty.Structures
using Magesty.Symmetries


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
input["regression"]["datafile"] = joinpath(@__DIR__, "EMBSET")

system = System(input, verbosity = false)

const NUM_CELLS = 27  # Center cell + 26 neighboring image cells
const NUM_FE = 32
const NUM_GE = 32
const NUM_ATOMS = NUM_FE + NUM_GE  # 64
const A_FEGE = 9.378              # Cubic lattice constant from input.toml

@testset "FeGe B20 2x2x2 Tests" begin
	@testset "Structure Tests" begin
		structure = system.structure

		# Lattice and reciprocal vectors
		expected_lattice = SMatrix{3, 3, Float64}([
			A_FEGE 0.0    0.0;
			0.0    A_FEGE 0.0;
			0.0    0.0    A_FEGE
		])
		@test structure.supercell.lattice_vectors ≈ expected_lattice atol = 1e-5

		expected_reciprocal = SMatrix{3, 3, Float64}([
			1/A_FEGE 0.0      0.0;
			0.0      1/A_FEGE 0.0;
			0.0      0.0      1/A_FEGE
		])
		@test structure.supercell.reciprocal_vectors ≈ expected_reciprocal atol = 1e-5

		# Basic properties
		@test structure.supercell.num_atoms == NUM_ATOMS
		@test structure.supercell.num_elements == 2  # Fe and Ge
		@test structure.supercell.kd_int_list ==
			  vcat(fill(1, NUM_FE), fill(2, NUM_GE))
		@test all(structure.is_periodic)

		# Atomic positions match input.toml
		expected_positions_vvec = input["structure"]["position"]
		expected_positions = reduce(hcat, expected_positions_vvec)
		@test structure.supercell.x_frac ≈ expected_positions atol = 1e-5

		@test structure.is_periodic == [true, true, true]
		@test structure.kd_name == ["Fe", "Ge"]
		@test structure.x_image_frac[:, :, 1] ≈ structure.supercell.x_frac atol = 1e-5
		@test structure.x_image_frac[:, :, 2] ≈ structure.supercell.x_frac .- 1.0 atol = 1e-5
		@test structure.x_image_frac[:, :, NUM_CELLS] ≈ structure.supercell.x_frac .+ 1.0 atol =
			1e-5
		@test structure.x_image_cart[:, :, 1] ≈
			  structure.supercell.lattice_vectors * structure.supercell.x_frac atol = 1e-5
		@test structure.x_image_cart[:, :, 2] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .- 1.0) atol =
			1e-5
		@test structure.x_image_cart[:, :, NUM_CELLS] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .+ 1.0) atol =
			1e-5
		@test structure.exist_image == [true for _ in 1:NUM_CELLS]
		# Atom-type grouping: first NUM_FE indices are Fe, next NUM_GE are Ge
		@test structure.atomtype_group == [collect(1:NUM_FE), collect((NUM_FE+1):NUM_ATOMS)]
	end

	@testset "Symmetry Tests" begin
		symmetry = system.symmetry
		structure = system.structure

		# Consistency checks (avoid hard-coding space-group specifics).
		@test symmetry.nat_prim ≥ 1
		@test symmetry.nsym ≥ 1
		@test symmetry.ntran ≥ 1
		@test symmetry.nat_prim * symmetry.ntran == NUM_ATOMS
		@test !isempty(symmetry.symdata)
		@test length(symmetry.atoms_in_prim) == symmetry.nat_prim

		# Identity operation must preserve atomic positions.
		identity_op = symmetry.symdata[1]
		for i in 1:structure.supercell.num_atoms
			pos = SVector{3, Float64}(structure.supercell.x_frac[:, i])
			transformed_pos = identity_op.rotation_frac * pos + identity_op.translation_frac
			@test transformed_pos ≈ pos atol = 1e-5
		end
	end

	sclus = SpinCluster(system, input, verbosity = false)
	Magesty.write_xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))

	@testset "SCE Regression Roundtrip (energy + torque)" begin
		spinconfigs = sclus.optimize.spinconfig_list
		# Use at most 2×(number of SALC groups) configs to keep test lightweight.
		num_cfg_target = 2 * length(sclus.basisset.salc_list)
		num_cfg = min(length(spinconfigs), num_cfg_target)
		spinconfigs = spinconfigs[1:num_cfg]

		num_atoms = sclus.structure.supercell.num_atoms
		salc_list = sclus.basisset.salc_list

		# Build design matrices (bias column included for energy, not for torque)
		design_E = Optimize.build_design_matrix_energy(salc_list, spinconfigs, sclus.symmetry)
		design_T = Optimize.build_design_matrix_torque(
			salc_list,
			spinconfigs,
			num_atoms,
			sclus.symmetry,
		)

		num_salcs = length(salc_list)
		@test size(design_E, 2) == num_salcs + 1
		@test size(design_T, 2) == num_salcs

		# Synthetic ground-truth coefficients and corresponding observations
		rng = Xoshiro(0x12345678)
		jphi_true = randn(rng, num_salcs)
		j0_true = 0.1234

		observed_energy_list = design_E[:, 2:end] * jphi_true .+ j0_true
		torque_flat = design_T * jphi_true
		block_size = 3 * num_atoms
		observed_torque_list = [
			reshape(torque_flat[((i-1)*block_size+1):(i*block_size)], 3, num_atoms) for
			i in 1:num_cfg
		]

		# Recover coefficients via elastic-net regression (torque-only weight).
		weight = 1.0
		j0_hat, jphi_hat = Optimize.elastic_net_regression(
			design_E,
			design_T,
			observed_energy_list,
			observed_torque_list,
			0.0,
			0.0,
			weight,
		)

		@test isapprox(j0_hat, j0_true; rtol = 1e-8, atol = 1e-8)
		@test isapprox(jphi_hat, jphi_true; rtol = 1e-8, atol = 1e-8)
	end
end
