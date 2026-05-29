using Test
using Random
using LinearAlgebra
using Statistics
using Magesty
using MultivariateStats: ridge as _mvs_ridge

# Component tests for the OLS / Ridge / AdaptiveRidge Cholesky path.
#
# The four properties exercised here are:
#
#   1. OLS Cholesky on a well-conditioned design agrees with an
#      independent pivoted-QR reference (`Matrix \\ Vector`) to ~1e-10.
#      Pins the algebraic identity `cholesky(X'X) \\ X'y == X \\ y`
#      modulo Cholesky-vs-QR rounding.
#
#   2. Rank-deficient design (duplicated column) triggers `PosDefException`
#      inside the OLS solver and is rethrown as `ArgumentError` whose
#      message names `Ridge` and explains the likely cause. This is the
#      explicit user signal that replaces the old pivoted-QR fallback.
#
#   3. Ridge Cholesky agrees with `MultivariateStats.ridge(...; bias = false)`
#      to ~1e-10. The cross-check pins the algebraic equivalence of
#      `cholesky(X'X + lambda*I) \\ X'y` and the externally maintained
#      reference implementation. `MultivariateStats` is a test-only
#      dependency, kept solely for this comparison.
#
#   4. `Ridge(lambda = 0.0)` delegates to the OLS Cholesky path: a
#      rank-deficient design therefore raises the same `ArgumentError`
#      rather than silently returning garbage.

# Build a well-conditioned `(X, y)` for a synthetic regression. n rows,
# p columns, p < n. Cond around 1e2 — comfortably inside Float64 precision
# for both Cholesky and QR.
function _well_conditioned_system(rng, n::Int, p::Int)
	X = randn(rng, n, p)
	beta_true = randn(rng, p)
	y = X * beta_true .+ 0.01 .* randn(rng, n)
	return X, y
end

# Build a rank-deficient `(X, y)` by duplicating the last column.
function _rank_deficient_system(rng, n::Int, p::Int)
	X = randn(rng, n, p)
	X[:, p] = X[:, p - 1]  # exact duplication
	beta_true = randn(rng, p)
	y = X * beta_true .+ 0.01 .* randn(rng, n)
	return X, y
end

@testset "_format_bytes (internal helper)" begin
	# Boundary cases for the byte-count pretty-printer.
	@test Magesty.Fitting._format_bytes(0) == "0 B"
	@test Magesty.Fitting._format_bytes(1023) == "1023 B"
	@test Magesty.Fitting._format_bytes(1024) == "1.0 KB"
	@test Magesty.Fitting._format_bytes(1536) == "1.5 KB"
	@test Magesty.Fitting._format_bytes(1024^2) == "1.0 MB"
	@test Magesty.Fitting._format_bytes(1024^3) == "1.0 GB"
	@test Magesty.Fitting._format_bytes(1024^4) == "1.0 TB"
	# TB fallback for very large inputs.
	@test Magesty.Fitting._format_bytes(5 * 1024^4) == "5.0 TB"
	# Negative input is an overflow signal upstream; assert rather than format.
	@test_throws ArgumentError Magesty.Fitting._format_bytes(-1)
end

@testset "OLS / Ridge Cholesky solver" begin

	@testset "OLS Cholesky ≈ pivoted-QR reference (well-conditioned)" begin
		rng = MersenneTwister(20260526)
		X, y = _well_conditioned_system(rng, 200, 12)
		b_chol = Magesty.Fitting.solve_coefficients(OLS(), X, y)
		b_qr   = X \ y
		@test maximum(abs.(b_chol .- b_qr)) < 1e-10
	end

	@testset "OLS PosDefException → ArgumentError pointing to Ridge" begin
		rng = MersenneTwister(20260526)
		X, y = _rank_deficient_system(rng, 100, 8)
		err = nothing
		try
			Magesty.Fitting.solve_coefficients(OLS(), X, y)
		catch e
			err = e
		end
		@test err isa ArgumentError
		# The message must name `Ridge` (the recommended remedy) and
		# describe the failure as linear dependence; the test pins both
		# tokens so a future message rewrite does not silently lose the
		# Ridge pointer.
		@test occursin("Ridge", err.msg)
		@test occursin("linearly dependent", err.msg)
	end

	@testset "Ridge Cholesky ≈ MultivariateStats.ridge" begin
		rng = MersenneTwister(20260526)
		X, y = _well_conditioned_system(rng, 200, 12)
		for lambda in (1e-6, 1e-3, 1.0, 10.0)
			b_chol = Magesty.Fitting.solve_coefficients(
				Ridge(lambda = lambda), X, y,
			)
			b_mvs = _mvs_ridge(X, y, lambda; bias = false)
			@test maximum(abs.(b_chol .- b_mvs)) < 1e-10
		end
	end

	@testset "Ridge constructor rejects negative / non-finite lambda" begin
		# The new Cholesky-Ridge path documents `lambda > 0` ⇒ SPD; a
		# negative or NaN `lambda` violates that and would surface as a
		# misleading OLS-style `ArgumentError`. The constructor now
		# rejects such inputs up front, matching `AdaptiveRidge`.
		@test_throws ArgumentError Ridge(lambda = -1e-6)
		@test_throws ArgumentError Ridge(lambda = NaN)
		@test_throws ArgumentError Ridge(lambda = -Inf)
		# Boundary: zero is allowed (delegates to OLS path).
		@test Ridge(lambda = 0.0) isa Ridge
		@test Ridge(lambda = 1e-6) isa Ridge
	end

	@testset "Ridge(λ ≈ 0) delegates to OLS Cholesky" begin
		rng = MersenneTwister(20260526)
		X, y = _rank_deficient_system(rng, 100, 8)
		err = nothing
		try
			Magesty.Fitting.solve_coefficients(Ridge(lambda = 0.0), X, y)
		catch e
			err = e
		end
		@test err isa ArgumentError
		@test occursin("Ridge", err.msg)
	end

end
