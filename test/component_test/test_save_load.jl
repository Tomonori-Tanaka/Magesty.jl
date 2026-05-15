using Test
using Magesty

# Fixtures -------------------------------------------------------------

const SL_DIMER_TOML = joinpath(@__DIR__, "..", "examples", "dimer", "input.toml")
const SL_FEPT_TOML =
	joinpath(@__DIR__, "..", "examples", "fept_tetragonal_2x2x2", "input.toml")
const SL_FEPT_EMBSET =
	joinpath(@__DIR__, "..", "examples", "fept_tetragonal_2x2x2", "EMBSET.dat")
const SL_BASELINE_DIR = joinpath(@__DIR__, "baselines")


# Structural + numerical equality of two SALC bases. After an XML
# round-trip the SALC basis is a fresh object (`===` does not hold), so
# the design-matrix `(l, m, site)` ordering is verified by comparing the
# salc_list contents element-by-element.
function _assert_salcbasis_equal(a, b)
	@test length(a.salc_list) == length(b.salc_list)
	for (group_a, group_b) in zip(a.salc_list, b.salc_list)
		@test length(group_a) == length(group_b)
		for (cbc_a, cbc_b) in zip(group_a, group_b)
			@test cbc_a.ls == cbc_b.ls
			@test cbc_a.Lf == cbc_b.Lf
			@test cbc_a.Lseq == cbc_b.Lseq
			@test cbc_a.atoms == cbc_b.atoms
			@test cbc_a.multiplicity == cbc_b.multiplicity
			@test isapprox(cbc_a.coefficient, cbc_b.coefficient, atol = 1e-10)
			@test size(cbc_a.coeff_tensor) == size(cbc_b.coeff_tensor)
			@test isapprox(cbc_a.coeff_tensor, cbc_b.coeff_tensor, atol = 1e-10)
		end
	end
	return nothing
end

# Tests ----------------------------------------------------------------

@testset "save / load" begin
	@testset "extension dispatch" begin
		basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
		@test_throws ArgumentError save(basis, "model.txt")
		@test_throws ArgumentError save(basis, "model")
		@test_throws ArgumentError load(SCEBasis, "basis.json")
		@test_throws ArgumentError load(SCEModel, "model.yaml")

		# An uppercase `.XML` extension is accepted.
		mktempdir() do dir
			path = joinpath(dir, "BASIS.XML")
			save(basis, path)
			@test isfile(path)
		end
	end

	@testset "SCEBasis round-trip" begin
		basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
		mktempdir() do dir
			path = joinpath(dir, "basis.xml")
			save(basis, path)
			reloaded = load(SCEBasis, path)

			@test reloaded isa SCEBasis
			@test reloaded.isotropy == basis.isotropy
			@test reloaded.symmetry.tol == basis.symmetry.tol
			@test reloaded.structure.supercell.num_atoms ==
				  basis.structure.supercell.num_atoms
			@test reloaded.structure.supercell.lattice_vectors ≈
				  basis.structure.supercell.lattice_vectors
			@test reloaded.structure.supercell.x_frac ≈ basis.structure.supercell.x_frac
			_assert_salcbasis_equal(basis.salcbasis, reloaded.salcbasis)
		end
	end

	@testset "SCEModel round-trip" begin
		basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
		n_salc = length(basis.salcbasis.salc_list)
		# Construct a model directly from synthetic coefficients — the
		# round-trip is about serialization, not about the fit.
		model = SCEModel(basis, -3.141592653589793, collect(1:n_salc) ./ 7)
		mktempdir() do dir
			path = joinpath(dir, "model.xml")
			save(model, path)
			reloaded = load(SCEModel, path)

			@test reloaded isa SCEModel
			@test reloaded.j0 ≈ model.j0
			@test reloaded.jphi ≈ model.jphi
			@test reloaded.basis.isotropy == basis.isotropy
			@test reloaded.basis.symmetry.tol == basis.symmetry.tol
			_assert_salcbasis_equal(model.basis.salcbasis, reloaded.basis.salcbasis)
		end
	end

	@testset "load SCEBasis from an SCEModel XML" begin
		basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
		n_salc = length(basis.salcbasis.salc_list)
		model = SCEModel(basis, 1.5, collect(1.0:n_salc))
		mktempdir() do dir
			model_path = joinpath(dir, "model.xml")
			save(model, model_path)

			# load(SCEBasis, ...) accepts an SCEModel XML, dropping <JPhi>.
			from_model = load(SCEBasis, model_path)
			@test from_model isa SCEBasis
			_assert_salcbasis_equal(basis.salcbasis, from_model.salcbasis)

			# load(SCEModel, ...) on a basis-only XML is a clear error.
			basis_path = joinpath(dir, "basis.xml")
			save(basis, basis_path)
			@test_throws ArgumentError load(SCEModel, basis_path)
		end
	end

	@testset "basis identity (design matrices) after reload" begin
		# The real round-trip guarantee: a dataset built from the reloaded
		# basis produces the same design matrices as the original, i.e.
		# the `(l, m, site)` ordering survives serialization. Equality is
		# `≈`, not exact: SALC coupling tensors are serialized at
		# `%.15e`, one digit short of full Float64 round-trip.
		basis = SCEBasis(SL_FEPT_TOML; verbosity = false)
		dataset = SCEDataset(basis, SL_FEPT_EMBSET)
		mktempdir() do dir
			path = joinpath(dir, "fept_basis.xml")
			save(basis, path)
			reloaded = load(SCEBasis, path)
			dataset_reloaded = SCEDataset(reloaded, SL_FEPT_EMBSET)

			@test dataset_reloaded.X_E ≈ dataset.X_E
			@test dataset_reloaded.X_T ≈ dataset.X_T
			# Observations come straight from the EMBSET file — exact.
			@test dataset_reloaded.y_E == dataset.y_E
			@test dataset_reloaded.y_T == dataset.y_T
		end
	end

	# Byte-diff regression against committed baselines: catches accidental
	# schema drift in the XML writer / reader. The check is a load -> re-save
	# round-trip — fully deterministic (no SALC eigensolve), so byte equality
	# holds regardless of platform.
	#
	# A "fresh build" byte-diff (`SCEBasis(toml)` -> save -> byte-compare to
	# committed baseline) is intentionally not done here: both fept and fege
	# have Lf=2 SALC groups that span 2D degenerate eigenvalue subspaces, and
	# LAPACK eigensolvers do not guarantee a platform-stable choice of basis
	# within such a subspace. Per-platform consistency of the build is still
	# verified by the "basis identity (design matrices) after reload" testset
	# above, which compares design matrices with `≈`.
	#
	# A principled fix (subspace-span comparison or canonical gauge fix at
	# SALC build) is tracked in `docs/design-notes/post-step7-cleanup.md`.
	@testset "byte-diff regression vs committed baseline" begin
		for (name, T) in (
			("fept_basis", SCEBasis),
			("fept_model", SCEModel),
			("fege_basis", SCEBasis),
		)
			baseline = joinpath(SL_BASELINE_DIR, "$name.xml")
			obj = load(T, baseline)
			mktempdir() do dir
				path = joinpath(dir, "$name.xml")
				save(obj, path)
				@test read(path, String) == read(baseline, String)
			end
		end
	end
end
