# Sunny.jl round-trip validation for the `sce_to_sunny` exporter.
#
# Sunny is a heavy dependency, so this lives in a dedicated environment
# (test/sunny/Project.toml) and runs only via `make test-sunny`; it is NOT part
# of `make test-all`. The check evaluates the generated script's system-building
# section in an actual Sunny session and confirms `energy(sys)` reproduces
# `predict_energy(model, sd) - model.j0` for random spin configurations — the
# correctness gate for the SCE → Sunny conversion (signs, normalization, the
# s = 1 single-ion factor, bond offsets, and primitive unfolding).

using Test
using Magesty
using Sunny
using Random
using LinearAlgebra: norm, inv
using StaticArrays: SVector

const _INTEG = joinpath(@__DIR__, "..", "integration")
_fix(parts...) = joinpath(_INTEG, parts...)

# Evaluate only the system-construction portion of a generated script (before the
# magnetic-structure block, which would randomize and minimize the state).
function _eval_system(script::AbstractString)
    head = split(script, "# ---- Magnetic structure")[1]
    m = Module(:GenSunnyScript)
    Core.eval(m, :(using Sunny))
    Core.eval(m, Meta.parseall(head))
    return Core.eval(m, :sys)
end

# Build the Sunny system from a generated script, map the model's supercell atoms
# onto Sunny sites (reshaping to the supercell for the primitive route), and
# return (is_primitive, max energy error over `n` random configs).
function _roundtrip(model::SCEModel; spin = 5 // 2, n::Int = 20)
    script = sce_to_sunny(model; spin = spin)
    primitive = occursin("Cell route: primitive", script)
    sys = _eval_system(script)

    L = Matrix{Float64}(model.basis.structure.supercell.lattice_vectors)
    xf = model.basis.structure.supercell.x_frac
    nat = size(xf, 2)

    if primitive
        prim = Magesty._sunny_build_primitive(model).prim
        target = reshape_supercell(sys, Matrix(prim.reshape_matrix))
        # Map each supercell atom to its Sunny site by matching the global
        # Cartesian position under supercell periodicity. (More robust at cell
        # boundaries than `position_to_site`, which is strict about ±1 wraps.)
        wrapcart(v) = (g = L \ v; L * (g .- round.(g)))
        ssites = collect(eachsite(target))
        spos = [Sunny.global_position(target, s) for s in ssites]
        sites = map(1:nat) do a
            ra = L * SVector{3, Float64}(xf[:, a])
            ssites[argmin([norm(wrapcart(spos[k] - ra)) for k in eachindex(ssites)])]
        end
    else
        target = sys
        sites = [(1, 1, 1, a) for a = 1:nat]
    end
    @test length(unique(sites)) == nat   # site mapping is a bijection

    rng = MersenneTwister(20260529)
    maxerr = 0.0
    for _ = 1:n
        sd = randn(rng, 3, nat)
        for c = 1:nat
            sd[:, c] ./= norm(sd[:, c])
        end
        for a = 1:nat
            set_dipole!(target, sd[:, a], sites[a])
        end
        maxerr = max(maxerr, abs(energy(target) - (Magesty.predict_energy(model, sd) - model.j0)))
    end
    return primitive, maxerr
end

@testset "Sunny round-trip" begin
    # Every fitted model unfolds onto the primitive cell. The 2×2×2 fixtures have
    # pairs with multiplicity > 1 (equal-distance degenerate bonds placed as
    # separate primitive bonds); fept also exercises the single-ion path. Tolerance
    # allows for the script's 12-significant-digit coupling literals.
    @testset "primitive route: $(f[1])" for f in
                                            (("dimer", "dimer.xml"),
        ("dimer", "dimer_dmi.xml"),
        ("chain", "chain.xml"),
        ("febcc_2x2x2_pm", "scecoeffs.xml"),
        ("fept_tetragonal_2x2x2", "scecoeffs.xml"),
        ("fege_2x2x2", "scecoeffs.xml"))
        model = Magesty.load(SCEModel, _fix(f...))
        primitive, err = _roundtrip(model)
        @test primitive
        @test err < 1e-8
    end

    # Energy invariance under the spin choice: the s_i s_j (bond) and the mode
    # factor (single-ion) cancel, so `energy(sys)` reproduces `predict_energy - j0`
    # for any (half-integer) spin. Sunny requires s to be a multiple of 1/2.
    @testset "energy invariant under spin: $(s)" for s in (1, 3 // 2, 5 // 2)
        model = Magesty.load(SCEModel, _fix("fept_tetragonal_2x2x2", "scecoeffs.xml"))
        _, err = _roundtrip(model; spin = s, n = 10)
        @test err < 1e-8
    end

    # Magnon dispersion scales as ħω ∝ 1/s (fixed energy landscape). From the same
    # aligned start both spins minimize to the same state, so the whole dispersion
    # of spin = 1 is 5/2 times that of spin = 5//2.
    @testset "dispersion scales as 1/s" begin
        model = Magesty.load(SCEModel, _fix("chain", "chain.xml"))
        function bandwidth(spin)
            sys = _eval_system(sce_to_sunny(model; spin = spin))
            for st in eachsite(sys)
                set_dipole!(sys, [0.0, 0.0, 1.0], st)
            end
            minimize_energy!(sys)
            swt = SpinWaveTheory(sys; measure = nothing)
            path = q_space_path(sys.crystal, [[0, 0, 0], [0.5, 0, 0], [0.25, 0.1, 0]], 16)
            return maximum(dispersion(swt, path))
        end
        @test bandwidth(1) ≈ 2.5 * bandwidth(5 // 2) rtol = 1e-5
    end

    # The explicit (folded supercell) route stays available on request and remains
    # exact.
    @testset "explicit route on request: $(f[1])" for f in
                                                      (("fept_tetragonal_2x2x2", "scecoeffs.xml"),)
        model = Magesty.load(SCEModel, _fix(f...))
        script = sce_to_sunny(model; spin = 5 // 2, placement = :explicit)
        @test occursin("Cell route: explicit", script)
        sys = _eval_system(script)
        nat = size(model.basis.structure.supercell.x_frac, 2)
        rng = MersenneTwister(7)
        maxerr = 0.0
        for _ = 1:10
            sd = randn(rng, 3, nat)
            for c = 1:nat
                sd[:, c] ./= norm(sd[:, c])
            end
            for a = 1:nat
                set_dipole!(sys, sd[:, a], (1, 1, 1, a))
            end
            maxerr = max(maxerr, abs(energy(sys) - (Magesty.predict_energy(model, sd) - model.j0)))
        end
        @test maxerr < 1e-8
    end
end
