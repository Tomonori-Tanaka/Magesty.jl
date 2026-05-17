using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
embset_path = joinpath(@__DIR__, "EMBSET")

basis = SCEBasis(input; verbosity = false)
spinconfigs_all = Magesty.read_embset(embset_path)

const NUM_CELLS = 27  # Center cell + 26 neighboring image cells
const NUM_FE = 32
const NUM_GE = 32
const NUM_ATOMS = NUM_FE + NUM_GE  # 64
const A_FEGE = 9.378              # Cubic lattice constant from input.toml

@testset "FeGe B20 2x2x2 Tests" begin
	@testset "Structure Tests" begin
		structure = basis.structure

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
		symmetry = basis.symmetry
		structure = basis.structure

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

	# Fit on a subset and persist the model — exercises the full new-API flow
	# on a large (64-atom) two-species system.
	num_cfg_target = 2 * length(basis.salcbasis.salc_list)
	num_cfg = min(length(spinconfigs_all), num_cfg_target)
	spinconfigs = spinconfigs_all[1:num_cfg]
	dataset = SCEDataset(basis, spinconfigs)
	fitted = fit(SCEFit, dataset, OLS(); torque_weight = 1.0, verbosity = false)
	Magesty.save(SCEModel(fitted), joinpath(@__DIR__, "scecoeffs.xml"))

	@testset "SCE Regression Roundtrip (energy + torque)" begin
		num_atoms = basis.structure.supercell.num_atoms
		num_salcs = length(basis.salcbasis.salc_list)

		@test size(dataset.X_E, 2) == num_salcs + 1
		@test size(dataset.X_T, 2) == num_salcs

		# Synthesize ground-truth coefficients and the matching observations.
		rng = Xoshiro(0x12345678)
		jphi_true = randn(rng, num_salcs)
		j0_true = 0.1234
		y_E_synth = dataset.X_E[:, 2:end] * jphi_true .+ j0_true
		y_T_synth = dataset.X_T * jphi_true

		synth_dataset = SCEDataset(
			basis,
			dataset.spinconfigs,
			dataset.X_E,
			dataset.X_T,
			y_E_synth,
			y_T_synth,
		)

		# Recover the coefficients via OLS (torque-only weight = 1.0).
		f = fit(SCEFit, synth_dataset, OLS(); torque_weight = 1.0, verbosity = false)
		@test isapprox(intercept(f), j0_true; rtol = 1e-8, atol = 1e-8)
		@test isapprox(coef(f), jphi_true; rtol = 1e-8, atol = 1e-8)
	end
end
