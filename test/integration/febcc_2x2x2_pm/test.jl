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

const NUM_CELLS = 27  # Total number of cells: center cell and its neighboring virtual cells
@testset "Fe BCC 2x2x2 Tests" begin
	@testset "Structure Tests" begin
		structure = basis.structure

		# Test lattice and basic properties
		a = 2.83  # Fe lattice constant in angstrom
		expected_lattice = SMatrix{3, 3, Float64}([
			2*a  0.0  0.0;
			0.0  2*a  0.0;
			0.0  0.0  2*a
		])
		@test structure.supercell.lattice_vectors ≈ expected_lattice atol=1e-5

		expected_reciprocal = SMatrix{3, 3, Float64}([
			1/(2*a)  0.0  0.0;
			0.0  1/(2*a)  0.0;
			0.0  0.0  1/(2*a)
		])
		@test structure.supercell.reciprocal_vectors ≈ expected_reciprocal atol=1e-5


		# Test basic properties
		@test structure.supercell.num_atoms == 16  # 2x2x2 supercell of BCC (2 atoms/unit cell)
		@test structure.supercell.num_elements == 1  # Only Fe atoms
		@test structure.supercell.kd_int_list == [1 for _ in 1:structure.supercell.num_atoms]
		@test all(structure.is_periodic)  # Periodic in all directions

		# Test atomic positions
		# BCC positions in 2x2x2 supercell
		expected_positions_vvec = [
			[0.00, 0.00, 0.00],
			[0.25, 0.25, 0.25],
			[0.00, 0.00, 0.50],
			[0.50, 0.00, 0.00],
			[0.00, 0.50, 0.00],
			[0.25, 0.25, 0.75],
			[0.75, 0.25, 0.25],
			[0.25, 0.75, 0.25],
			[0.00, 0.50, 0.50],
			[0.50, 0.00, 0.50],
			[0.50, 0.50, 0.00],
			[0.25, 0.75, 0.75],
			[0.75, 0.25, 0.75],
			[0.75, 0.75, 0.25],
			[0.50, 0.50, 0.50],
			[0.75, 0.75, 0.75],
		]
		expected_positions = reduce(hcat, expected_positions_vvec)
		@test structure.supercell.x_frac ≈ expected_positions atol=1e-5

		@test structure.is_periodic == [true, true, true]
		@test structure.kd_name == ["Fe"]
		@test structure.x_image_frac[:, :, 1] ≈ structure.supercell.x_frac atol=1e-5
		@test structure.x_image_frac[:, :, 2] ≈ structure.supercell.x_frac .- 1.0 atol=1e-5
		@test structure.x_image_frac[:, :, NUM_CELLS] ≈ structure.supercell.x_frac .+ 1.0 atol=1e-5
		@test structure.x_image_cart[:, :, 1] ≈
			  structure.supercell.lattice_vectors * structure.supercell.x_frac atol=1e-5
		@test structure.x_image_cart[:, :, 2] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .- 1.0) atol=1e-5
		@test structure.x_image_cart[:, :, NUM_CELLS] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .+ 1.0) atol=1e-5
		@test structure.exist_image == [true for _ in 1:NUM_CELLS]
		@test structure.atomtype_group ==
			  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]
	end

	@testset "Symmetry Tests" begin
		symmetry = basis.symmetry

		# Test basic properties
		@test symmetry.international_symbol == "Im-3m"  # Space group for BCC Fe
		@test symmetry.spacegroup_number == 229
		@test symmetry.nat_prim == 1  # Two atoms in primitive cell
		@test symmetry.nsym == 768  # 48 point group operations × 16 translations for 2x2x2
		@test symmetry.ntran == 16   # Number of translations for 2x2x2 supercell
		@test !isempty(symmetry.symdata)
		@test length(symmetry.atoms_in_prim) == symmetry.nat_prim

		# Test that identity operation preserves atomic positions
		structure = basis.structure
		identity_op = symmetry.symdata[1]
		for i in 1:structure.supercell.num_atoms
			pos = SVector{3, Float64}(structure.supercell.x_frac[:, i])
			transformed_pos = identity_op.rotation_frac * pos + identity_op.translation_frac
			@test transformed_pos ≈ pos atol=1e-5
		end
	end

	# Fit on a subset and persist the model to XML — exercises the full
	# SCEBasis → SCEDataset → SCEFit → save round-trip on a realistic system.
	num_cfg_target = 2 * length(basis.salcbasis.salc_list)
	num_cfg = min(length(spinconfigs_all), num_cfg_target)
	spinconfigs = spinconfigs_all[1:num_cfg]
	dataset = SCEDataset(basis, spinconfigs)
	fitted = fit(SCEFit, dataset, OLS(); torque_weight = 0.5, verbosity = false)
	Magesty.save(SCEModel(fitted), joinpath(@__DIR__, "scecoeffs.xml"))

	@testset "SCE Regression Roundtrip (energy + torque)" begin
		num_atoms = basis.structure.supercell.num_atoms
		num_salcs = length(basis.salcbasis.salc_list)

		@test size(dataset.X_E, 2) == num_salcs
		@test size(dataset.X_T, 2) == num_salcs

		# Synthesize ground-truth coefficients and the matching observations.
		rng = Xoshiro(0x12345678)
		jphi_true = randn(rng, num_salcs)
		j0_true = 0.1234
		y_E_synth = dataset.X_E * jphi_true .+ j0_true
		y_T_synth = dataset.X_T * jphi_true

		synth_dataset = SCEDataset(
			basis,
			dataset.spinconfigs,
			dataset.X_E,
			dataset.X_T,
			y_E_synth,
			y_T_synth,
		)

		# Recover the coefficients across the full torque-weight range.
		for w in [0.0, 0.5, 1.0]
			f = fit(SCEFit, synth_dataset, OLS(); torque_weight = w, verbosity = false)
			@test isapprox(intercept(f), j0_true; rtol = 1e-8, atol = 1e-8)
			@test isapprox(coef(f), jphi_true; rtol = 1e-8, atol = 1e-8)
		end
	end

	@testset "AdaptiveLasso smoke + XML round-trip" begin
		# AdaptiveLasso runs end-to-end on the same realistic SCE
		# pipeline (basis -> dataset -> fit -> save -> load). The OLS
		# pilot default is fine on this fixture because the design is
		# well-conditioned at this configuration count. A small lambda
		# is chosen on purpose: this is a numerical-safety / round-trip
		# smoke test, not a support-recovery test.
		f_ada = fit(
			SCEFit, dataset, AdaptiveLasso(lambda = 1e-4);
			torque_weight = 0.5, verbosity = false,
		)
		@test all(isfinite, coef(f_ada))
		@test isfinite(intercept(f_ada))

		mktempdir() do dir
			xml_path = joinpath(dir, "scecoeffs_adalasso.xml")
			model = SCEModel(f_ada)
			Magesty.save(model, xml_path)
			reloaded = Magesty.load(SCEModel, xml_path)
			@test reloaded.j0 ≈ model.j0
			@test reloaded.jphi ≈ model.jphi
		end
	end

	@testset "AdaptiveRidge smoke + XML round-trip" begin
		# AdaptiveRidge runs the iterative L0-approximation end-to-end on
		# the same realistic SCE pipeline (basis -> dataset -> fit ->
		# save -> load). A small lambda is chosen on purpose: this is a
		# numerical-safety / round-trip smoke test, not a support-recovery
		# test.
		f_ar = fit(
			SCEFit, dataset, AdaptiveRidge(lambda = 1e-4);
			torque_weight = 0.5, verbosity = false,
		)
		@test all(isfinite, coef(f_ar))
		@test isfinite(intercept(f_ar))

		mktempdir() do dir
			xml_path = joinpath(dir, "scecoeffs_adaridge.xml")
			model = SCEModel(f_ar)
			Magesty.save(model, xml_path)
			reloaded = Magesty.load(SCEModel, xml_path)
			@test reloaded.j0 ≈ model.j0
			@test reloaded.jphi ≈ model.jphi
		end
	end

	@testset "AdaptiveLasso(::SCEFit) / AdaptiveLasso(::SCEModel) reuse" begin
		# Convenience constructors must wrap the prior fit's coefficients
		# in a PrecomputedPilot and forward the rest of the kwargs to the
		# underlying keyword ctor. `fitted` is the OLS SCEFit defined
		# above; it shares the same dataset and SCEBasis as the
		# AdaptiveLasso call below, so the length check inside
		# `solve_coefficients(::PrecomputedPilot, X, y)` is satisfied.
		est_from_fit = AdaptiveLasso(fitted; lambda = 1e-4)
		@test est_from_fit isa AdaptiveLasso
		@test est_from_fit.pilot isa PrecomputedPilot
		@test est_from_fit.pilot.beta == coef(fitted)
		@test est_from_fit.lambda == 1e-4
		@test est_from_fit.gamma == 1.0  # default

		model_ols = SCEModel(fitted)
		est_from_model = AdaptiveLasso(model_ols; lambda = 1e-4, gamma = 0.5)
		@test est_from_model.pilot isa PrecomputedPilot
		@test est_from_model.pilot.beta == coef(model_ols)
		@test est_from_model.gamma == 0.5

		# End-to-end dispatch: pilot regression is skipped, the stored
		# coefficients drive the adaptive weights, and the resulting
		# coefficients are finite.
		f_reuse = fit(
			SCEFit, dataset, est_from_fit;
			torque_weight = 0.5, verbosity = false,
		)
		@test all(isfinite, coef(f_reuse))
		@test isfinite(intercept(f_reuse))

		# Passing `pilot = ...` on top of the positional SCEFit would
		# otherwise be silently overridden by Julia's kwarg-splat
		# semantics. The convenience constructor guards against that
		# footgun explicitly with an ArgumentError.
		@test_throws ArgumentError AdaptiveLasso(
			fitted; pilot = OLS(), lambda = 1e-4,
		)
	end
end
