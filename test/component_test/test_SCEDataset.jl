using Test
using Magesty
using AtomsBase
using Unitful

# Fixtures -------------------------------------------------------------

const DIMER_TOML_DS = joinpath(@__DIR__, "..", "examples", "dimer", "input.toml")

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
        @test size(d.X_E) == (2, n_salc + 1)
        @test all(d.X_E[:, 1] .== 1.0)            # bias column
        @test size(d.X_T) == (2 * block, n_salc)  # no bias column
        @test d.y_E == [-1.0, 1.0]
        @test d.y_T == vcat(vec(configs[1].torques), vec(configs[2].torques))
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
        d1.X_E[1, 1] = -999.0
        @test d.X_E[1, 1] == 1.0

        # reordering reorders both X_E rows and X_T blocks consistently
        dr = d[[2, 1]]
        @test dr.X_E == d.X_E[[2, 1], :]
        @test dr.X_T == d.X_T[vcat((block + 1):(2 * block), 1:block), :]

        # boolean indexing
        db = d[[true, false]]
        @test length(db) == 1
        @test db.X_E == d.X_E[1:1, :]
    end

    @testset "vcat reconstructs and checks basis identity" begin
        d = SCEDataset(basis, configs)
        merged = vcat(d[1], d[2])
        @test length(merged) == 2
        @test merged.X_E == d.X_E
        @test merged.X_T == d.X_T
        @test merged.y_E == d.y_E
        @test merged.y_T == d.y_T
        @test merged.basis === basis

        # datasets built from different SCEBasis objects cannot be vcat-ed
        other_basis = SCEBasis(DIMER_TOML_DS; verbosity = false)
        d_other = SCEDataset(other_basis, configs)
        @test_throws ArgumentError vcat(d, d_other)
    end
end
