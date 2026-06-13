using Test
using Magesty
using AtomsBase
using Unitful

# Fixtures -------------------------------------------------------------

const DIMER_TOML_DS = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

# FM / AFM spin configurations for the 2-atom dimer.
function _dimer_configs()
    SC = Magesty.SpinConfigs.SpinConfig
    fm = SC(-1.0, [1.0, 1.0],
            [0.0 0.0; 0.0 0.0; 1.0 1.0],
            [0.0 0.0; 0.0 0.0; 0.0 0.0])
    afm = SC(1.0, [1.0, 1.0],
             [0.0 0.0; 0.0 0.0; -1.0 1.0],
             [0.0 0.0; 0.0 0.0; 0.0 0.0])
    return [fm, afm]
end

# A minimal periodic 2-atom system for the AtomsBase sugar constructor.
function _periodic_pair_ds()
    pos = ([0.0, 0.0, 0.0]u"Å", [1.5, 1.5, 1.5]u"Å")
    atoms = [Atom(:Fe, p) for p in pos]
    cv = ([3.0, 0.0, 0.0]u"Å", [0.0, 3.0, 0.0]u"Å", [0.0, 0.0, 3.0]u"Å")
    return periodic_system(atoms, cv)
end

# Tests ----------------------------------------------------------------

@testset "SCEDataset" begin
    basis = SCEBasis(DIMER_TOML_DS; verbosity = false)
    configs = _dimer_configs()
    n_atoms = basis.structure.supercell.num_atoms
    block = 3 * n_atoms

    @testset "base constructor builds design matrices" begin
        d = SCEDataset(basis, configs)
        @test d isa SCEDataset
        @test d.basis === basis
        @test length(d) == 2
        @test lastindex(d) == 2
        n_salc = length(basis.salcbasis.salc_list)
        @test size(d.X_E) == (2, n_salc)          # no bias column
        @test size(d.X_T) == (2 * block, n_salc)  # no bias column
        @test d.y_E == [-1.0, 1.0]
        @test d.y_T == vcat(vec(configs[1].torques), vec(configs[2].torques))
    end

    @testset "rejects zero moment on a basis-referenced atom" begin
        # The dimer's 2-body SALC references both atoms. A configuration that
        # puts a zero moment on a referenced atom is ill-posed: its undefined
        # spin direction would feed the energy basis and bias the fit, so the
        # dataset constructor rejects it.
        SC = Magesty.SpinConfigs.SpinConfig
        sc = SC(-1.0, [0.0, 1.0],
                [0.0 0.0; 0.0 0.0; 1.0 1.0],
                [0.0 0.0; 0.0 0.0; 0.0 0.0])
        @test_throws ArgumentError SCEDataset(basis, [sc])
    end

    @testset "embset_path constructor matches explicit configs" begin
        # write a temp EMBSET mirroring the FM / AFM configs above
        embset = tempname()
        open(embset, "w") do io
            println(io, "# config 1 (FM)")
            println(io, "-1.0")
            println(io, "1  0.0 0.0  1.0  0.0 0.0 0.0")
            println(io, "2  0.0 0.0  1.0  0.0 0.0 0.0")
            println(io, "# config 2 (AFM)")
            println(io, "1.0")
            println(io, "1  0.0 0.0 -1.0  0.0 0.0 0.0")
            println(io, "2  0.0 0.0  1.0  0.0 0.0 0.0")
        end
        d_file = SCEDataset(basis, embset)
        d_ref = SCEDataset(basis, configs)
        @test length(d_file) == 2
        @test d_file.X_E ≈ d_ref.X_E
        @test d_file.X_T ≈ d_ref.X_T
        @test d_file.y_E == d_ref.y_E
        rm(embset)
    end

    @testset "TOML sugar constructor" begin
        d = SCEDataset(DIMER_TOML_DS, configs; verbosity = false)
        d_ref = SCEDataset(basis, configs)
        @test d isa SCEDataset
        @test size(d.X_E) == size(d_ref.X_E)
        @test size(d.X_T) == size(d_ref.X_T)
        # a fresh basis is built internally — not the same object
        @test d.basis !== basis
    end

    @testset "AtomsBase sugar constructor" begin
        sys = _periodic_pair_ds()
        d = SCEDataset(sys, configs;
            interaction = (body1 = (lmax = Dict(:Fe => 2),),
                           body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
        @test d isa SCEDataset
        @test length(d) == 2
        @test size(d.X_E, 1) == 2
        @test size(d.X_T, 1) == 2 * block
    end

    @testset "getindex copies and slices rows / blocks" begin
        d = SCEDataset(basis, configs)

        d1 = d[1]
        @test d1 isa SCEDataset
        @test length(d1) == 1
        @test d1.X_E == d.X_E[1:1, :]
        @test d1.X_T == d.X_T[1:block, :]
        @test d1.y_E == d.y_E[1:1]

        d2 = d[2:2]
        @test d2.X_E == d.X_E[2:2, :]
        @test d2.X_T == d.X_T[(block + 1):(2 * block), :]   # second config's block

        # getindex returns a copy: mutating the slice leaves the parent intact
        original_first_entry = d.X_E[1, 1]
        d1.X_E[1, 1] = -999.0
        @test d.X_E[1, 1] == original_first_entry

        # reordering reorders both X_E rows and X_T blocks consistently
        dr = d[[2, 1]]
        @test dr.X_E == d.X_E[[2, 1], :]
        @test dr.X_T == d.X_T[vcat((block + 1):(2 * block), 1:block), :]

        # boolean indexing
        db = d[[true, false]]
        @test length(db) == 1
        @test db.X_E == d.X_E[1:1, :]
    end

    @testset "SCEModel / SCEFit reuse the embedded basis" begin
        d_ref = SCEDataset(basis, configs)
        # `Ridge(lambda = 1e-8)` rather than `OLS()` because the dimer
        # FM / AFM configs have spins along z, which makes every torque
        # zero by S × S = 0; at the default `torque_weight = 1.0` the
        # design matrix is then entirely zero. The test only checks
        # basis reuse, so the regularizer is incidental.
        f = fit(SCEFit, d_ref, Ridge(lambda = 1e-8); verbosity = false)
        model = SCEModel(f)

        d_from_model = SCEDataset(model, configs)
        @test d_from_model isa SCEDataset
        @test d_from_model.basis === basis   # basis reused, not rebuilt
        @test d_from_model.X_E == d_ref.X_E
        @test d_from_model.X_T == d_ref.X_T

        d_from_fit = SCEDataset(f, configs)
        @test d_from_fit isa SCEDataset
        @test d_from_fit.basis === basis
        @test d_from_fit.X_E == d_ref.X_E
        @test d_from_fit.X_T == d_ref.X_T
    end

    @testset "vcat checks SALC-basis fingerprint compatibility" begin
        d = SCEDataset(basis, configs)
        merged = vcat(d[1], d[2])
        @test length(merged) == 2
        @test merged.X_E == d.X_E
        @test merged.X_T == d.X_T
        @test merged.y_E == d.y_E
        @test merged.y_T == d.y_T
        @test merged.basis === basis

        # A basis rebuilt from the same input is a distinct object but
        # has an identical structural fingerprint (the fingerprint is a
        # deterministic hash of the SALC structure). Its dataset can be
        # vcat-ed: the design-matrix column orderings provably agree.
        twin_basis = SCEBasis(DIMER_TOML_DS; verbosity = false)
        @test twin_basis !== basis
        @test twin_basis.salc_fingerprint == basis.salc_fingerprint
        d_twin = SCEDataset(twin_basis, configs)
        merged_twin = vcat(d, d_twin)
        @test length(merged_twin) == 4
        @test merged_twin.basis === basis              # keeps the first basis
        @test merged_twin.X_E == vcat(d.X_E, d_twin.X_E)
        @test merged_twin.X_T == vcat(d.X_T, d_twin.X_T)

        # A basis with a different SALC structure has a different
        # fingerprint, so its dataset cannot be vcat-ed.
        d_diff = SCEDataset(_periodic_pair_ds(), configs;
            interaction = (body1 = (lmax = Dict(:Fe => 2),),
                           body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
        @test d_diff.basis.salc_fingerprint != basis.salc_fingerprint
        @test_throws ArgumentError vcat(d, d_diff)
    end
end
