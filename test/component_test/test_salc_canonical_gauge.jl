using Test
using LinearAlgebra
using Magesty

const _canonicalize_eigenspace = Magesty.SALCBases._canonicalize_eigenspace

# `_canonicalize_eigenspace(V)` is the gauge-fix used inside SALC build to
# collapse the LAPACK-implementation-dependent choice of eigenbasis within
# a degenerate eigenspace. The contract:
#   - the output `W` has orthonormal columns,
#   - `W` spans the same subspace as `V` (projectors agree),
#   - any orthonormal basis of the same subspace produces the same `W`
#     up to round-off (axis-ordered modified Gram-Schmidt on `P * eⱼ`),
#   - signs are pinned by "first nonzero axis projection is positive".

@testset "SALCBases _canonicalize_eigenspace" begin
    @testset "orthonormality and subspace preservation (2D in R^4)" begin
        V = [1.0  0.0;
             0.0  1.0;
             1.0  0.0;
             0.0  1.0] ./ sqrt(2)
        W = _canonicalize_eigenspace(V)
        @test size(W) == size(V)
        @test isapprox(W' * W, I(2); atol = 1e-12)
        @test isapprox(W * W', V * V'; atol = 1e-12)
    end

    @testset "gauge invariance under right-multiplication by O(d)" begin
        V = [1.0  0.0;
             0.0  1.0;
             1.0  0.0;
             0.0  1.0] ./ sqrt(2)
        # Rotate within the column space by an arbitrary angle.
        θ = 0.37
        Q = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        W1 = _canonicalize_eigenspace(V)
        W2 = _canonicalize_eigenspace(V * Q)
        @test isapprox(W1, W2; atol = 1e-12)

        # A reflection (det Q = -1) also leaves the subspace invariant.
        R = [1.0 0.0; 0.0 -1.0]
        W3 = _canonicalize_eigenspace(V * R)
        @test isapprox(W1, W3; atol = 1e-12)
    end

    @testset "axis-aligned subspace yields the axis basis" begin
        # S = span(e1, e2) in R^4, with V rotated within S by 30 degrees.
        θ = π / 6
        V = [cos(θ) -sin(θ);
             sin(θ)  cos(θ);
             0.0     0.0;
             0.0     0.0]
        W = _canonicalize_eigenspace(V)
        @test isapprox(W, [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]; atol = 1e-12)
    end

    @testset "axis-order pivoting skips axes orthogonal to S" begin
        # S = span(e2, e3); e1 and e4 are orthogonal to S, so P*e1 and
        # P*e4 vanish and the canonical basis falls back to (e2, e3).
        V = [0.0       0.0;
             1.0       1.0;
             1.0      -1.0;
             0.0       0.0] ./ sqrt(2)
        W = _canonicalize_eigenspace(V)
        expected = [0.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 0.0]
        @test isapprox(W, expected; atol = 1e-12)
    end

    @testset "1D subspace: sign pinned by first nonzero axis component" begin
        V_pos = reshape([3.0 / 5, 4.0 / 5, 0.0], 3, 1)
        @test isapprox(_canonicalize_eigenspace(V_pos),
                       reshape([3.0 / 5, 4.0 / 5, 0.0], 3, 1);
                       atol = 1e-12)

        # Sign-flipped input describes the same subspace, so the output
        # must be identical to the positive case.
        V_neg = reshape([-3.0 / 5, -4.0 / 5, 0.0], 3, 1)
        @test isapprox(_canonicalize_eigenspace(V_neg),
                       reshape([3.0 / 5, 4.0 / 5, 0.0], 3, 1);
                       atol = 1e-12)
    end

    @testset "tol failure raises when no axis projection is large enough" begin
        V = reshape([0.6, 0.8, 0.0], 3, 1)
        # Maximum axis projection norm is 0.8 < tol = 1.5, so no column
        # survives Gram-Schmidt and the helper must error.
        @test_throws ErrorException _canonicalize_eigenspace(V; tol = 1.5)
    end
end
