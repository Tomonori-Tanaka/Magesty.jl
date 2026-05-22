using Test
using Magesty

# Direct access to the symmetry-reduction primitives. `Cluster` is exported
# from `Magesty`, but `is_translationally_equiv_cluster` and the
# `cluster_dict` / `irreducible_cluster_dict` / `cluster_orbits_dict`
# fields are internal contracts — qualified access makes that explicit.
const _Clusters = Magesty.Clusters
const _Structures = Magesty.Structures
const _Symmetries = Magesty.Symmetries
const _InputSpecs = Magesty.InputSpecs

# Build (structure, symmetry, cluster) directly from a TOML-style dict.
# Cluster is built inside SCEBasis but discarded; for unit tests we
# reproduce the same call path without going through SALC construction.
function _build_test_cluster(input::AbstractDict)
	system, interaction, options = _InputSpecs.parse_toml_inputs(input)
	structure = _Structures.Structure(system; verbosity = false)
	symmetry = _Symmetries.Symmetry(structure, options; verbosity = false)
	cluster = _Clusters.Cluster(structure, symmetry, interaction; verbosity = false)
	return structure, symmetry, cluster
end

# Irreducible representatives for `body`, as sorted Vector{Int} clusters.
# `SortedCounter` iterates over keys in sorted order; the keys here are
# already sorted atom-index lists from `irreducible_clusters`.
function _irreducible_reps(cluster, body::Integer)
	return Vector{Int}[copy(c) for c in cluster.irreducible_cluster_dict[body]]
end

# Flatten the per-body orbit dict into a single list of clusters.
function _orbit_clusters(cluster, body::Integer)
	flat = Vector{Vector{Int}}()
	for (_, orbit) in cluster.cluster_orbits_dict[body]
		append!(flat, orbit)
	end
	return flat
end

# Structural invariants that hold for every lattice and every body order.
# Pinning these — instead of pinning exact orbit counts — lets the test
# catch broken symmetry reductions without depending on supercell-wrap
# arithmetic that may itself be non-obvious.
function _check_invariants(symmetry, cluster, body::Integer)
	reps = _irreducible_reps(cluster, body)
	orbits = cluster.cluster_orbits_dict[body]
	orbit_flat = _orbit_clusters(cluster, body)

	# (1) Self-equivalence: every cluster is translationally equivalent to itself.
	for c in reps
		@test _Clusters.is_translationally_equiv_cluster(c, c, symmetry)
	end

	# (2) Symmetry of the equivalence relation on distinct representatives.
	for i in eachindex(reps), j in (i+1):lastindex(reps)
		c1, c2 = reps[i], reps[j]
		@test _Clusters.is_translationally_equiv_cluster(c1, c2, symmetry) ==
			  _Clusters.is_translationally_equiv_cluster(c2, c1, symmetry)
	end

	# (3) Distinct irreducible representatives must NOT be translationally
	# equivalent — that is exactly the property `irreducible_clusters`
	# is supposed to enforce.
	for i in eachindex(reps), j in (i+1):lastindex(reps)
		@test !_Clusters.is_translationally_equiv_cluster(reps[i], reps[j], symmetry)
	end

	# (4) Orbits partition the irreducible representatives:
	#     union = full set, no cluster appears in two orbits.
	@test Set(orbit_flat) == Set(reps)
	@test length(orbit_flat) == length(reps)  # no duplicates

	return nothing
end

# ----------------------------------------------------------------------
# Lattice fixtures
# ----------------------------------------------------------------------

# Simple cubic, 1 atom in the primitive cell, 2x2x2 supercell.
# Lattice constant a = 1.0 (primitive) → supercell side = 2.0.
# 1NN distance = a = 1.0, 2NN distance = sqrt(2) ≈ 1.414. Cutoff 1.2
# isolates 1NN only.
#
# Analytical prediction (Pm-3m / Oh):
# - 1NN coordination = 6 directions (±x, ±y, ±z), but in a 2x2x2 cell
#   ±dir wrap to the same supercell atom, so each prim atom has 3
#   unique NN targets.
# - Pure-translation classes within cluster_dict[2][prim]: 3 (one per
#   axis direction; translation cannot rotate).
# - Full space-group orbit: 1 (Oh maps x↔y↔z).
const _SC_INPUT = Dict(
	"general" => Dict(
		"name" => "sc",
		"kd" => ["X"],
		"nat" => 8,
		"periodicity" => [true, true, true],
	),
	"symmetry" => Dict("tolerance" => 1e-5),
	"interaction" => Dict(
		"nbody" => 2,
		"body1" => Dict("lmax" => Dict("X" => 0)),
		"body2" => Dict("cutoff" => Dict("X-X" => 1.2), "lsum" => 2),
	),
	"structure" => Dict(
		"kd_list" => fill(1, 8),
		"lattice" => [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
		"position" => [
			[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
			[0.0, 0.5, 0.0], [0.5, 0.5, 0.0],
			[0.0, 0.0, 0.5], [0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5], [0.5, 0.5, 0.5],
		],
	),
)

# BCC conventional cell, 2-atom basis, 2x2x2 supercell.
# Conventional lattice constant = 5.66, primitive = 1 atom in BCC.
# NN distance = a_conv * sqrt(3)/2 ≈ 2.451; 2NN = a_conv ≈ 2.83.
# Cutoff 2.5 isolates 1NN only.
#
# Analytical prediction (Im-3m / Oh + body-centering):
# - 1NN coordination = 8 body-center directions per atom.
# - Pure-translation classes within cluster_dict[2][prim]: 4. The body-
#   centering translation (1/2, 1/2, 1/2) plus the inverse mapping makes
#   each <±±±>/2 direction equivalent to its antipode; 8 directions / 2
#   = 4 translation classes.
# - Full space-group orbit: 1 (Oh acts transitively on <111>/2).
const _BCC_INPUT = Dict(
	"general" => Dict(
		"name" => "bcc",
		"kd" => ["Y"],
		"nat" => 16,
		"periodicity" => [true, true, true],
	),
	"symmetry" => Dict("tolerance" => 1e-5),
	"interaction" => Dict(
		"nbody" => 2,
		"body1" => Dict("lmax" => Dict("Y" => 0)),
		"body2" => Dict("cutoff" => Dict("Y-Y" => 2.5), "lsum" => 2),
	),
	"structure" => Dict(
		"kd_list" => fill(1, 16),
		"lattice" => [[5.66, 0.0, 0.0], [0.0, 5.66, 0.0], [0.0, 0.0, 5.66]],
		"position" => [
			[0.00, 0.00, 0.00], [0.25, 0.25, 0.25],
			[0.00, 0.00, 0.50], [0.50, 0.00, 0.00],
			[0.00, 0.50, 0.00], [0.25, 0.25, 0.75],
			[0.75, 0.25, 0.25], [0.25, 0.75, 0.25],
			[0.00, 0.50, 0.50], [0.50, 0.00, 0.50],
			[0.50, 0.50, 0.00], [0.25, 0.75, 0.75],
			[0.75, 0.25, 0.75], [0.75, 0.75, 0.25],
			[0.50, 0.50, 0.50], [0.75, 0.75, 0.75],
		],
	),
)

# 2D square lattice, 1-atom primitive, 4x4x1 supercell with a 10× vacuum
# along z (no periodicity in z). Lattice side along x and y = 4.0.
# 1NN distance = 1.0; diagonal NN = sqrt(2). Cutoff 1.2 isolates 1NN.
#
# Analytical prediction (p4mm in plane, C4v point group):
# - 1NN coordination = 4 in-plane directions (±x, ±y).
# - Pure-translation classes within cluster_dict[2][prim]: 2. Within a
#   4×4 supercell, ±x and ±y are distinct supercell atoms; but the
#   translation that sends +x_neighbor → prim_atom also sends prim_atom
#   → -x_neighbor, so {prim, +x} ~ {prim, -x} under translation
#   (similarly for ±y). x-axis and y-axis classes stay distinct under
#   pure translation.
# - Full space-group orbit: 1 (C4 rotates x↔y).
const _SQUARE_INPUT = Dict(
	"general" => Dict(
		"name" => "square",
		"kd" => ["Z"],
		"nat" => 16,
		"periodicity" => [true, true, false],
	),
	"symmetry" => Dict("tolerance" => 1e-5, "isotropy" => true),
	"interaction" => Dict(
		"nbody" => 2,
		"body1" => Dict("lmax" => Dict("Z" => 0)),
		"body2" => Dict("cutoff" => Dict("Z-Z" => 1.2), "lsum" => 2),
	),
	"structure" => Dict(
		"kd_list" => fill(1, 16),
		"lattice" => [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 10.0]],
		"position" => [
			[0.00, 0.00, 0.0], [0.25, 0.00, 0.0],
			[0.50, 0.00, 0.0], [0.75, 0.00, 0.0],
			[0.00, 0.25, 0.0], [0.25, 0.25, 0.0],
			[0.50, 0.25, 0.0], [0.75, 0.25, 0.0],
			[0.00, 0.50, 0.0], [0.25, 0.50, 0.0],
			[0.50, 0.50, 0.0], [0.75, 0.50, 0.0],
			[0.00, 0.75, 0.0], [0.25, 0.75, 0.0],
			[0.50, 0.75, 0.0], [0.75, 0.75, 0.0],
		],
	),
)

# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

@testset "_translation_canonical_form" begin
	@testset "Simple cubic 2x2x2" begin
		_, symmetry, cluster = _build_test_cluster(_SC_INPUT)
		reps = _irreducible_reps(cluster, 2)
		canonical = _Clusters._translation_canonical_form

		# (1) Sort-invariance: shuffling the input does not change the
		# canonical form.
		for c in reps
			@test canonical(c, symmetry) == canonical(reverse(c), symmetry)
		end

		# (2) Orbit consistency: every translation (and inverse) of c
		# shares the same canonical form as c.
		for c in reps
			base = canonical(c, symmetry)
			for sym_tran in symmetry.symnum_translation
				translated = [symmetry.map_sym[atom, sym_tran] for atom in c]
				@test canonical(translated, symmetry) == base
				translated_inv = [symmetry.map_sym_inv[atom, sym_tran] for atom in c]
				@test canonical(translated_inv, symmetry) == base
			end
		end

		# (3) Distinct irreducible representatives carry distinct
		# canonical forms.
		canonical_forms = [canonical(r, symmetry) for r in reps]
		@test length(unique(canonical_forms)) == length(reps)

		# (4) Cross-check against the predicate: equiv iff same canonical.
		for i in eachindex(reps), j in eachindex(reps)
			equiv = _Clusters.is_translationally_equiv_cluster(reps[i], reps[j], symmetry)
			same_canonical = canonical(reps[i], symmetry) == canonical(reps[j], symmetry)
			@test equiv == same_canonical
		end
	end

	@testset "BCC conventional 2x2x2" begin
		_, symmetry, cluster = _build_test_cluster(_BCC_INPUT)
		reps = _irreducible_reps(cluster, 2)
		canonical = _Clusters._translation_canonical_form

		canonical_forms = [canonical(r, symmetry) for r in reps]
		@test length(unique(canonical_forms)) == length(reps)

		for c in reps
			base = canonical(c, symmetry)
			for sym_tran in symmetry.symnum_translation
				translated = [symmetry.map_sym[atom, sym_tran] for atom in c]
				@test canonical(translated, symmetry) == base
			end
		end
	end

	@testset "Square 4x4" begin
		_, symmetry, cluster = _build_test_cluster(_SQUARE_INPUT)
		reps = _irreducible_reps(cluster, 2)
		canonical = _Clusters._translation_canonical_form

		canonical_forms = [canonical(r, symmetry) for r in reps]
		@test length(unique(canonical_forms)) == length(reps)
	end
end

@testset "Clusters: symmetry reduction" begin
	@testset "Simple cubic 2x2x2, 1NN only" begin
		_, symmetry, cluster = _build_test_cluster(_SC_INPUT)
		_check_invariants(symmetry, cluster, 2)

		# Pinned (Pm-3m): in a 2x2x2 supercell ±dir wrap to the same atom,
		# so each primitive atom has 3 unique NN targets (x/y/z axes).
		# All 3 are space-group-equivalent under Oh -> a single orbit
		# whose size equals the irreducible count.
		@test length(cluster.irreducible_cluster_dict[2]) == 3
		@test length(cluster.cluster_orbits_dict[2]) == 1
		@test length(cluster.cluster_orbits_dict[2][1]) == 3
	end

	@testset "BCC conventional 2x2x2, 1NN only" begin
		_, symmetry, cluster = _build_test_cluster(_BCC_INPUT)
		_check_invariants(symmetry, cluster, 2)

		# Pinned (Im-3m): 8 <111>/2 NN directions per atom; the body-
		# centering pure translation identifies each direction with its
		# antipode, leaving 4 translation-irreducible representatives.
		# Oh acts transitively on <111>/2 -> a single space-group orbit.
		@test length(cluster.irreducible_cluster_dict[2]) == 4
		@test length(cluster.cluster_orbits_dict[2]) == 1
		@test length(cluster.cluster_orbits_dict[2][1]) == 4
	end

	@testset "Square 4x4, 1NN only" begin
		_, symmetry, cluster = _build_test_cluster(_SQUARE_INPUT)
		_check_invariants(symmetry, cluster, 2)

		# Pinned (p4mm in plane): 4 in-plane NN directions per atom; in a
		# 4x4 supercell ±dir wrap to distinct atoms, but translation
		# identifies +x with -x and +y with -y, leaving 2 translation-
		# irreducible representatives. C4 maps x↔y -> a single orbit.
		@test length(cluster.irreducible_cluster_dict[2]) == 2
		@test length(cluster.cluster_orbits_dict[2]) == 1
		@test length(cluster.cluster_orbits_dict[2][1]) == 2
	end
end
