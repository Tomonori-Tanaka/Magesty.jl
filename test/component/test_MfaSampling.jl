"""
Unit tests for the code-agnostic `Magesty.MfaSampling` sampler and the
`sample_mfa_incar` VASP orchestration.

The von Mises-Fisher draw is validated against its analytic mean resultant
length ``A_3(κ) = \\coth κ - 1/κ`` rather than pinned to a recorded sample, so
the test certifies the hand-rolled p=3 sampler reproduces the distribution.
Deterministic behavior (magnitude preservation, fixed atoms, the low-τ ordered
limit) is checked exactly.
"""

using Magesty
using Test
using LinearAlgebra
using Random
using StaticArrays

using Magesty.MfaSampling:
    thermal_averaged_m, tau_from_magnetization, sample_vmf_direction,
    mfa_sample, mfa_sweep, parse_atom_index_spec

@testset "von Mises-Fisher mean resultant length A_3(κ)" begin
    Random.seed!(20260605)
    mean_dir = SVector(0.0, 0.0, 1.0)
    N = 200_000
    for κ in (0.5, 1.0, 5.0, 20.0)
        acc = SVector(0.0, 0.0, 0.0)
        for _ = 1:N
            acc += sample_vmf_direction(mean_dir, κ)
        end
        mean_vec = acc / N
        Rbar = norm(mean_vec)
        A3 = coth(κ) - 1 / κ
        @test isapprox(Rbar, A3; atol = 5e-3)
        # Mean direction aligns with mean_dir (transverse components vanish).
        @test isapprox(mean_vec[1], 0.0; atol = 5e-3)
        @test isapprox(mean_vec[2], 0.0; atol = 5e-3)
        @test mean_vec[3] > 0
    end
    # Sampled directions are unit vectors.
    for κ in (0.0, 1e-6, 3.0)
        d = sample_vmf_direction(mean_dir, κ)
        @test isapprox(norm(d), 1.0; atol = 1e-12)
    end
end

@testset "MFA self-consistency round trip" begin
    # Boundary limits (analytic): the Langevin self-consistency m = L(3m/τ)
    # gives m -> 1 as τ -> 0 (fully ordered) and m -> 0 as τ -> 1 (disordered).
    @test thermal_averaged_m(1e-7) == 1.0
    @test thermal_averaged_m(2.0) == 0.0
    @test tau_from_magnetization(0.0) == 1.0
    @test tau_from_magnetization(1.0) == 0.0
    # m(τ) is monotonically decreasing in τ on (0, 1).
    @test thermal_averaged_m(0.3) > thermal_averaged_m(0.7)
    # Inverse consistency: τ -> m -> τ recovers the input.
    for τ in (0.2, 0.5, 0.8)
        m = thermal_averaged_m(τ)
        @test isapprox(tau_from_magnetization(m), τ; atol = 1e-6)
    end
    # Small / large magnetizations whose τ lies outside the bracket are clamped
    # rather than crashing the bracketed solver.
    @test tau_from_magnetization(1e-4) ≥ 0.99
    @test tau_from_magnetization(1 - 1e-4) ≤ 0.01
end

@testset "mfa_sample magnitude preservation and limits" begin
    Random.seed!(1)
    # Columns: atom1 = [0,0,2], atom2 = [3,0,0], atom3 = [0,0,0] (zero-norm).
    spins = [0.0 3.0 0.0; 0.0 0.0 0.0; 2.0 0.0 0.0]
    mags = [norm(spins[:, i]) for i = 1:3]

    out = mfa_sample(spins, "tau", 0.5)
    for i = 1:3
        @test isapprox(norm(out[:, i]), mags[i]; atol = 1e-10)
    end
    # Zero-norm column stays zero.
    @test all(out[:, 3] .== 0.0)

    # τ below MIN_TEMP returns the input unchanged (fully ordered limit).
    @test mfa_sample(spins, "tau", 1e-7) == spins

    # "m" path maps onto the same machinery and preserves magnitudes.
    out_m = mfa_sample(spins, "m", 0.6)
    for i = 1:3
        @test isapprox(norm(out_m[:, i]), mags[i]; atol = 1e-10)
    end

    @test_throws ArgumentError mfa_sample(spins, "bogus", 0.5)
end

@testset "mfa_sweep count, fixed, and uniform atoms" begin
    Random.seed!(2)
    spins = [0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0;
             2.0 2.0 2.0 2.0]  # 4 atoms, all +z, magnitude 2
    configs = mfa_sweep(spins; variable = "tau", start = 0.1, stop = 0.5,
        num_points = 3, num_samples = 2)
    @test length(configs) == 3 * 2
    @test all(size(c) == (3, 4) for c in configs)

    # Fixed atoms (no randomize) are byte-identical to the input column.
    cfgs_fixed = mfa_sweep(spins; variable = "tau", start = 0.4, stop = 0.4,
        num_points = 1, num_samples = 4, fixed_indices = [1, 3])
    for c in cfgs_fixed
        @test c[:, 1] == spins[:, 1]
        @test c[:, 3] == spins[:, 3]
    end

    # Uniform atoms keep magnitude but are (almost surely) reoriented.
    cfgs_uniform = mfa_sweep(spins; variable = "tau", start = 0.05, stop = 0.05,
        num_points = 1, num_samples = 1, uniform_indices = [2])
    @test isapprox(norm(cfgs_uniform[1][:, 2]), 2.0; atol = 1e-10)

    # Validation.
    @test_throws ArgumentError mfa_sweep(spins; variable = "m", start = 0.5,
        stop = 0.4, num_points = 2)  # start > stop
    @test_throws ArgumentError mfa_sweep(spins; variable = "tau", start = 0.1,
        stop = 0.5, num_points = 0)
    @test_throws ArgumentError mfa_sweep(spins; variable = "tau", start = 0.1,
        stop = 0.5, num_points = 2, fixed_indices = [5])  # out of range
end

@testset "randomize draws a Haar-uniform rotation" begin
    # Exact identities of the Haar measure on SO(3), used as the oracle:
    #   <R> = 0            (invariance under left multiplication)
    #   <tr R> = 0         (the rotation angle satisfies <cos θ> = -1/2 under
    #                       p(θ) = (1 - cos θ)/π, and tr R = 1 + 2 cos θ)
    # Per-entry variance is <R_ij^2> = 1/3, so σ(<R_ij>) = 1/sqrt(3N); for
    # tr R, Var(cos θ) = 1/4 gives σ(<tr R>) = 1/sqrt(N).
    #
    # The pre-fix sampler returned the minimal rotation carrying +z to a random
    # axis. That family has <R_xx> = <R_yy> = 1/2 and hence <tr R> = 1, so both
    # gates below fail on it by >100σ — they resolve the defect with margin, not
    # merely a fitted tolerance.
    Random.seed!(11)
    N = 20_000
    Rsum = zeros(3, 3)
    for _ = 1:N
        Rsum .+= Magesty.MfaSampling._random_rotation_matrix()
    end
    Rmean = Rsum / N
    @test isapprox(Rmean, zeros(3, 3); atol = 0.03)   # 7σ headroom (σ ≈ 0.0041)
    @test isapprox(tr(Rmean), 0.0; atol = 0.05)       # 7σ headroom (σ ≈ 0.0071)

    # Every sampled rotation must be a proper rotation.
    Random.seed!(12)
    for _ = 1:200
        R = Magesty.MfaSampling._random_rotation_matrix()
        @test isapprox(R' * R, I; atol = 1e-12)
        @test isapprox(det(R), 1.0; atol = 1e-12)
    end
end

@testset "randomize is isotropic regardless of the reference orientation" begin
    # <n n^T> = I/3 is exact for any isotropic distribution of unit vectors, and
    # constrains the polar and azimuthal marginals at once. τ below MIN_TEMP puts
    # the sweep in the ordered branch, so the vMF scatter is switched off and only
    # the global rotation is under test.
    #
    # The reference is deliberately non-collinear and not aligned with +z: the
    # pre-fix sampler was isotropic for a +z reference by construction, and biased
    # the azimuth by ~8x for an in-plane one (deviation 0.17 in <n n^T> for the x
    # column). Diagonal entries have σ = sqrt(4/45)/sqrt(N) ≈ 0.0021 at N = 20000,
    # so atol = 0.02 leaves ~9σ of headroom while the defect exceeds it 8-fold.
    Random.seed!(13)
    N = 20_000
    refs = [SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0),
        normalize(SVector(1.0, 1.0, 1.0))]
    spins = reduce(hcat, refs)
    configs = mfa_sweep(spins; variable = "tau", start = 1.0e-6, stop = 1.0e-6,
        num_points = 1, num_samples = N, randomize = true)
    @test length(configs) == N
    for atom = 1:length(refs)
        C = zeros(3, 3)
        for c in configs
            n = SVector{3, Float64}(c[1, atom], c[2, atom], c[3, atom])
            C .+= n * n'
        end
        C ./= N
        @test isapprox(C, Matrix(I / 3, 3, 3); atol = 0.02)
    end
end

@testset "parse_atom_index_spec" begin
    @test parse_atom_index_spec(""; max_index = 10) == Int[]
    @test parse_atom_index_spec("1-10,12,20-22"; max_index = 22) ==
          vcat(1:10, 12, 20:22)
    @test parse_atom_index_spec("3 1 5-7"; max_index = 10) == [1, 3, 5, 6, 7]
    @test parse_atom_index_spec("4:2"; max_index = 10) == [2, 3, 4]  # reversed range
    @test_throws ArgumentError parse_atom_index_spec("1-x"; max_index = 10)
    @test_throws ArgumentError parse_atom_index_spec("0"; max_index = 10)
    @test_throws ArgumentError parse_atom_index_spec("11"; max_index = 10)
end

@testset "sample_mfa_incar end to end" begin
    Random.seed!(3)
    incar = joinpath(@__DIR__, "fixtures", "incar", "INCAR")
    mktempdir() do dir
        outdir = joinpath(dir, "out")
        paths = sample_mfa_incar(incar; variable = "tau", start = 0.1, stop = 0.3,
            num_points = 3, num_samples = 2, outdir = outdir, fix = "1")
        # num_points * num_samples files, numbered 1..6.
        @test length(paths) == 6
        @test basename.(paths) == ["sample-$(i).INCAR" for i = 1:6]
        @test all(isfile, paths)

        input = Magesty.IncarIO.parse_incar(incar)
        atom1_in = input[:MAGMOM][1:3]
        for p in paths
            d = Magesty.IncarIO.parse_incar(p)
            m = d[:MAGMOM]
            # MAGMOM and M_CONSTR are written identically.
            @test d[:M_CONSTR] == m
            # Per-atom magnitudes are preserved (input atoms all have |m| = 3).
            for i = 1:4
                @test isapprox(norm(m[3i-2:3i]), 3.0; atol = 1e-6)
            end
            # Fixed atom 1 is unchanged from the input.
            @test isapprox(m[1:3], atom1_in; atol = 1e-6)
            # Non-sampled keys are preserved.
            @test d[:ENCUT] == 520.0
        end
    end

    @test_throws ArgumentError sample_mfa_incar(incar; variable = "bogus",
        start = 0.1, stop = 0.3, num_points = 3)
end

println("All MfaSampling tests passed.")
