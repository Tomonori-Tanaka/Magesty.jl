using Test
using LinearAlgebra
using TOML
using Magesty

# Direct access to the projection-side primitives. `Cluster` is exported
# through `Magesty.Clusters`; everything else here is internal and
# reached through qualified module names.
const _SALC = Magesty.SALCBases
const _CB = Magesty.CoupledBases
const _IS_ZERO = _SALC.is_obviously_zero_coupled_basis_product
const _PROJECTION = _SALC.projection_matrix_coupled_basis

# Build a CoupledBasis whose metadata matches what the predicate looks
# at. The coefficient tensor is filled with zeros because
# `is_obviously_zero_coupled_basis_product` only inspects `(Lf, ls,
# atoms)`; the tensor never enters the cheap check.
function _make_cb(; ls::Vector{Int}, Lf::Int, atoms::Vector{Int})
	N = length(ls)
	@assert N == length(atoms)
	Lseq = N > 2 ? zeros(Int, N - 2) : Int[]
	site_dims = ntuple(i -> 2 * ls[i] + 1, N)
	tensor_dims = (site_dims..., 2 * Lf + 1)
	tensor = zeros(Float64, tensor_dims)
	return _CB.CoupledBasis(ls, Lf, Lseq, atoms, tensor)
end

@testset "SALCBases projection layer" begin
	# ---- is_obviously_zero_coupled_basis_product: pure metadata predicate

	@testset "is_obviously_zero_coupled_basis_product" begin
		base = _make_cb(ls = [1, 1], Lf = 0, atoms = [1, 2])

		# Matching metadata: the predicate is conservative — the actual
		# product may still be zero, but the cheap predicate must not
		# declare it zero up front.
		same = _make_cb(ls = [1, 1], Lf = 0, atoms = [1, 2])
		@test _IS_ZERO(base, same) == false
		@test _IS_ZERO(same, base) == false  # symmetric

		# A change in any of the three metadata fields short-circuits to
		# `true`. Each is independent, so test each in isolation.
		diff_Lf = _make_cb(ls = [1, 1], Lf = 2, atoms = [1, 2])
		@test _IS_ZERO(base, diff_Lf) == true
		@test _IS_ZERO(diff_Lf, base) == true

		diff_ls = _make_cb(ls = [1, 2], Lf = 0, atoms = [1, 2])
		@test _IS_ZERO(base, diff_ls) == true
		@test _IS_ZERO(diff_ls, base) == true

		diff_atoms = _make_cb(ls = [1, 1], Lf = 0, atoms = [1, 3])
		@test _IS_ZERO(base, diff_atoms) == true
		@test _IS_ZERO(diff_atoms, base) == true
	end

	# ---- projection_matrix_coupled_basis: invariants on a real fixture

	@testset "projection_matrix invariants (dimer fixture)" begin
		# Re-run the pipeline up to coupled-basis classification so we can
		# call `projection_matrix_coupled_basis` directly on each group.
		toml = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")
		input = TOML.parsefile(toml)
		system, interaction, options = Magesty.InputSpecs.parse_toml_inputs(input)
		structure = Magesty.Structures.Structure(system; verbosity = false)
		symmetry = Magesty.Symmetries.Symmetry(structure, options; verbosity = false)
		cluster = Magesty.Clusters.Cluster(
			structure, symmetry, interaction; verbosity = false,
		)

		classified = _SALC.construct_and_classify_coupled_basislist(
			structure, symmetry, cluster,
			interaction.body1_lmax, interaction.bodyn_lsum, interaction.nbody;
			isotropy = options.isotropy,
		)
		classified = _SALC.filter_basisdict(classified, symmetry)

		# Sanity: the dimer test fixture produces at least one non-empty
		# group; otherwise the assertions below would all be vacuous.
		@test !isempty(classified)

		for (key, group) in classified
			P = _PROJECTION(group, symmetry)
			d = size(P, 1)

			# (1) Symmetric. The matrix is built as the average of real
			# representation matrices, then accumulated to itself in the
			# inner loop. After the loop, it must satisfy P = P^T at
			# machine precision.
			@test isapprox(P, P'; atol = 1e-12)

			# (2) Idempotent: P^2 = P. The fundamental invariant of a
			# projector — `_compute_salc_groups` relies on it when it
			# rounds the eigenvalues at `digits = 6` and asserts they
			# are 0 or 1.
			@test isapprox(P * P, P; atol = 1e-10)

			# (3) Tolerance sweep on the idempotence error. Each
			# tolerance is satisfied by Float64 round-off in a properly
			# normalised projector; failing the tighter one signals an
			# accumulating numerical error in the averaging step.
			err = opnorm(P * P - P)
			@test err < 1e-6
			@test err < 1e-8
			@test err < 1e-10

			# (4) Eigenvalues lie at 0 or 1 (the only spectrum a
			# projector can have). `1e-6` is loose enough to absorb the
			# `digits = 6` rounding `_compute_salc_groups` applies
			# downstream.
			eigvals_P = real.(eigen(Symmetric(P)).values)
			for ev in eigvals_P
				@test isapprox(ev, 0.0; atol = 1e-6) ||
					  isapprox(ev, 1.0; atol = 1e-6)
			end

			# (5) Trace equals the rank (number of unit eigenvalues),
			# rounded to the nearest integer. Trace deviation pinpoints
			# a normalisation slip even when the eigenvalues round
			# cleanly to {0, 1}.
			rank_from_eigs = count(ev -> isapprox(ev, 1.0; atol = 1e-6), eigvals_P)
			@test isapprox(tr(P), rank_from_eigs; atol = 1e-8)
			@test 0 <= rank_from_eigs <= d
		end
	end
end
