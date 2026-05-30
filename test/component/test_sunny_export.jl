# Component tests for the SCE → Sunny.jl exporter decomposition (Sunny-free).
#
# The decomposition is the correctness core: it must reproduce
# `predict_energy(model, sd) - model.j0` exactly from the per-bond 3×3 exchange
# matrices and per-site single-ion matrices. The tesseral-harmonic conventions
# baked into the helper matrices are cross-checked against
# `TesseralHarmonics.Zₗₘ`.

using Test
using Magesty
using Random
using LinearAlgebra: norm, dot, tr, transpose
using StaticArrays: SVector

const _Z = Magesty.TesseralHarmonics.Zₗₘ
const _CBC = Magesty.CoupledBases.CoupledBasis_with_coefficient
_fixture(parts...) = joinpath(@__DIR__, "..", "integration", parts...)

_rand_dir(rng) = (v = randn(rng, 3); SVector{3, Float64}(v ./ norm(v)))

@testset "SunnyExport" begin
	rng = MersenneTwister(20260529)

	@testset "l=1 pair matrix matches tesseral contraction" begin
		# folded index 1,2,3 ↔ m = -1,0,+1.
		for _ = 1:50
			folded = randn(rng, 3, 3)
			M = Magesty._sunny_l1_pair_matrix(folded)
			ea = _rand_dir(rng)
			eb = _rand_dir(rng)
			ref = 0.0
			for m1 = -1:1, m2 = -1:1
				ref += folded[m1 + 2, m2 + 2] * _Z(1, m1, ea) * _Z(1, m2, eb)
			end
			@test dot(ea, M * eb) ≈ ref atol = 1e-12
		end
	end

	@testset "l=2 onsite matrix matches tesseral contraction" begin
		# folded index 1..5 ↔ m = -2..+2.
		for _ = 1:50
			folded = randn(rng, 5)
			A = Magesty._sunny_l2_onsite_matrix(folded)
			@test A ≈ transpose(A) atol = 1e-12          # symmetric
			@test tr(A) ≈ 0.0 atol = 1e-12               # traceless
			e = _rand_dir(rng)
			ref = 0.0
			for m = -2:2
				ref += folded[m + 3] * _Z(2, m, e)
			end
			@test dot(e, A * e) ≈ ref atol = 1e-12
		end
	end

	# Reconstruction equals predict_energy - j0 for any model whose only
	# spin-dependent SALCs are l=l=1 pairs and l=2 single-ion (all current
	# fixtures qualify; assert no skipped terms before comparing).
	function _check_reconstruction(path; n = 30, atol = 1e-9)
		model = Magesty.load(Magesty.SCEModel, path)
		terms = Magesty._sunny_build_terms(model)
		@test isempty(terms.skipped)
		nat = size(model.basis.structure.supercell.x_frac, 2)
		maxerr = 0.0
		for _ = 1:n
			sd = Matrix{Float64}(undef, 3, nat)
			for c = 1:nat
				sd[:, c] = _rand_dir(rng)
			end
			e_rec = Magesty._sunny_reconstruct_energy(terms, sd)
			e_pred = Magesty.predict_energy(model, sd) - model.j0
			maxerr = max(maxerr, abs(e_rec - e_pred))
		end
		@test maxerr < atol
		return terms
	end

	@testset "dimer Heisenberg reconstruction" begin
		terms = _check_reconstruction(_fixture("dimer", "dimer.xml"))
		@test length(terms.pairs) == 1
		@test isempty(terms.onsites)
		M = terms.pairs[(1, 2)]
		# Pure isotropic exchange: M ≈ J·I.
		@test M[1, 2] ≈ 0 atol = 1e-10
		@test M[1, 3] ≈ 0 atol = 1e-10
		@test M[2, 3] ≈ 0 atol = 1e-10
		@test M[1, 1] ≈ M[2, 2] atol = 1e-10
		@test M[2, 2] ≈ M[3, 3] atol = 1e-10
	end

	@testset "dimer DM reconstruction" begin
		terms = _check_reconstruction(_fixture("dimer", "dimer_dmi.xml"))
		M = terms.pairs[(1, 2)]
		# Pure DM: the symmetric part vanishes (M is antisymmetric).
		@test M + transpose(M) ≈ zeros(3, 3) atol = 1e-10
	end

	@testset "fept single-ion + anisotropic reconstruction" begin
		terms = _check_reconstruction(
			_fixture("fept_tetragonal_2x2x2", "scecoeffs.xml"); atol = 1e-7)
		@test !isempty(terms.onsites)
		for (_, A) in terms.onsites
			@test A ≈ transpose(A) atol = 1e-10
			@test tr(A) ≈ 0.0 atol = 1e-10
		end
	end

	@testset "unsupported SALCs are skipped" begin
		# Synthetic key groups: a valid l=l=1 pair, a 3-body cluster, and an
		# l=2 pair. Only the first is exportable.
		valid_ct = zeros(Float64, 3, 3, 1)
		valid_ct[3, 3, 1] = 1.0
		valid = _CBC([1, 1], 0, Int[], [1, 2], valid_ct, [1.0], 1, [[1, 2]])

		three = _CBC([1, 1, 1], 0, [1], [1, 2, 3],
			zeros(Float64, 3, 3, 3, 1), [1.0], 1, [[1, 2, 3]])

		l2pair = _CBC([2, 2], 0, Int[], [1, 2],
			zeros(Float64, 5, 5, 1), [1.0], 1, [[1, 2]])

		# Match the abstract element typing of a real `salc_list`
		# (`Vector{Vector{CoupledBasis_with_coefficient}}`).
		salc_list = Vector{_CBC}[_CBC[valid], _CBC[three], _CBC[l2pair]]
		terms = Magesty._sunny_decompose(salc_list, [1.0, 1.0, 1.0])

		@test haskey(terms.pairs, (1, 2))
		@test isempty(terms.onsites)
		@test length(terms.skipped) == 2
		@test any(s -> occursin("n_body=3", s), terms.skipped)
		@test any(s -> occursin("ls=[2, 2]", s), terms.skipped)
	end

	@testset "primitive build unfolds every fixture" begin
		# Magesty keeps only minimum-distance images (distant pairs = 0), so every
		# fitted model maps cleanly onto the primitive cell — including the minimal
		# 2×2×2 fixtures whose pairs have multiplicity > 1 (equal-distance
		# degenerate bonds, placed as separate primitive bonds).
		for f in (("dimer", "dimer.xml"), ("dimer", "dimer_dmi.xml"), ("chain", "chain.xml"),
			("febcc_2x2x2_pm", "scecoeffs.xml"), ("fept_tetragonal_2x2x2", "scecoeffs.xml"),
			("fege_2x2x2", "scecoeffs.xml"))
			pm = Magesty._sunny_build_primitive(Magesty.load(Magesty.SCEModel, _fixture(f...)))
			@test pm.clean
			@test !isempty(pm.bonds)
			@test isempty(pm.skipped)
		end
	end

	@testset "sce_to_sunny script generation" begin
		m = Magesty.load(Magesty.SCEModel, _fixture("dimer", "dimer_dmi.xml"))
		s = sce_to_sunny(m)
		@test Meta.parseall(s) isa Expr                       # syntactically valid
		@test occursin("using Sunny", s)
		@test occursin("Cell route: primitive", s)
		@test occursin("set_exchange!", s)
		@test occursin("SpinWaveTheory", s)
		@test occursin("q_space_path", s)
		# Plotting snippet: per-band lines, not the broken matrix-as-vertices call.
		@test !occursin("lines(disp')", s)
		@test occursin("for b in axes(disp, 1)", s)
		# q-point labels on the x-axis, energy unit on the y-axis, and thin vlines
		# at the interior high-symmetry path corners.
		@test occursin("xticks = path.xticks", s)
		@test occursin("ylabel = \"Energy (eV)\"", s)
		@test occursin("vlines!(ax, path.xticks[1][2:end-1]", s)

		# fept (multiplicity > 1, with single-ion) now auto-selects the unfolded
		# primitive route.
		me = Magesty.load(Magesty.SCEModel, _fixture("fept_tetragonal_2x2x2", "scecoeffs.xml"))
		sp = sce_to_sunny(me)
		@test occursin("Cell route: primitive", sp)
		@test occursin("set_onsite_coupling!(sys, S ->", sp)   # single-ion path
		@test Meta.parseall(sp) isa Expr

		# The explicit (folded supercell) route remains available on request.
		se = sce_to_sunny(me; placement = :explicit)
		@test occursin("Cell route: explicit", se)
		@test occursin("set_exchange_at!", se)
		@test occursin("set_onsite_coupling_at!", se)
		@test Meta.parseall(se) isa Expr

		# File output and argument validation.
		mktempdir() do dir
			p = joinpath(dir, "lswt")
			sce_to_sunny(m; output = p)
			@test isfile(p * ".jl")
		end
		@test_throws ArgumentError sce_to_sunny(m; placement = :nonsense)
	end
end
