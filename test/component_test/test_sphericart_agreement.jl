using Test
using LinearAlgebra
using StaticArrays
using SpheriCart
using ..MySphericalHarmonics
using ..MySphericalHarmonics: Zₗₘ_unsafe, ∂ᵢZlm_unsafe

# SpheriCart index: lm2idx(l,m) = m + l + l^2 + 1  (1-based)
_lm2idx(l, m) = m + l + l * l + 1

@testset "SpheriCart vs MySphericalHarmonics agreement" begin
    lmax = 4
    sph  = SphericalHarmonics(lmax)

    rng_pts = [
        normalize([1.0,  0.0,  0.0]),
        normalize([0.0,  1.0,  0.0]),
        normalize([0.0,  0.0,  1.0]),
        normalize([1.0,  2.0,  3.0]),
        normalize([-1.0, 1.0, -2.0]),
        normalize([0.3, -0.7,  0.5]),
        normalize([0.0,  1.0,  1.0]),
    ]

    @testset "Zₗₘ values" begin
        for uvec in rng_pts
            uvec_sv = SVector{3,Float64}(uvec...)
            Y_sph   = compute(sph, uvec_sv)   # SVector{(lmax+1)^2}

            for l in 0:lmax, m in -l:l
                y_my  = Zₗₘ_unsafe(l, m, uvec)
                y_sph = Y_sph[_lm2idx(l, m)]
                @test isapprox(y_my, y_sph; atol = 1e-12) broken=false
            end
        end
    end

    @testset "∂ᵢZlm gradients" begin
        for uvec in rng_pts
            uvec_sv = SVector{3,Float64}(uvec...)
            _, ∇Y_sph = compute_with_gradients(sph, uvec_sv)

            for l in 0:lmax, m in -l:l
                g_my  = ∂ᵢZlm_unsafe(l, m, uvec)
                g_sph = ∇Y_sph[_lm2idx(l, m)]
                @test isapprox(g_my[1], g_sph[1]; atol = 1e-10)
                @test isapprox(g_my[2], g_sph[2]; atol = 1e-10)
                @test isapprox(g_my[3], g_sph[3]; atol = 1e-10)
            end
        end
    end
end
