using Test
using Magesty
using AtomsBase
using Unitful
import TOML

# Fixtures -------------------------------------------------------------

const DIMER_TOML = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

# A minimal periodic 2-atom system, used for the AtomsBase / sublabel /
# fan-out / odd-warning tests.
function _periodic_pair(; atom_names = (nothing, nothing))
    pos = ([0.0, 0.0, 0.0]u"Å", [1.5, 1.5, 1.5]u"Å")
    atoms = map(zip(pos, atom_names)) do (p, nm)
        nm === nothing ? Atom(:Fe, p) : Atom(:Fe, p; atom_name = nm)
    end
    cv = ([3.0, 0.0, 0.0]u"Å", [0.0, 3.0, 0.0]u"Å", [0.0, 0.0, 3.0]u"Å")
    return periodic_system(collect(atoms), cv)
end

# Tests ----------------------------------------------------------------

@testset "SCEBasis" begin
    @testset "dict and TOML paths agree" begin
        d = TOML.parsefile(DIMER_TOML)
        b_dict = SCEBasis(d; verbosity = false)
        b_toml = SCEBasis(DIMER_TOML; verbosity = false)
        @test b_dict isa SCEBasis
        @test b_toml isa SCEBasis
        @test length(b_dict.salcbasis.salc_list) == length(b_toml.salcbasis.salc_list)
        @test b_dict.structure.supercell.num_atoms == b_toml.structure.supercell.num_atoms
        @test b_dict.structure.kd_name == b_toml.structure.kd_name
    end

    @testset "keyword path matches TOML path" begin
        b_toml = SCEBasis(DIMER_TOML; verbosity = false)
        b_kw = SCEBasis(;
            lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
            kd = [:X], kd_list = [1, 1],
            positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1]],
            periodicity = (false, false, false),
            interaction = (body1 = (lmax = Dict(:X => 0),),
                           body2 = (lsum = 2, cutoff = Dict((:X, :X) => -1.0))),
            isotropy = true, verbosity = false,
        )
        @test length(b_kw.salcbasis.salc_list) == length(b_toml.salcbasis.salc_list)
        @test b_kw.structure.supercell.num_atoms == b_toml.structure.supercell.num_atoms
        @test b_kw.structure.kd_name == b_toml.structure.kd_name
        @test b_kw.structure.supercell.x_frac ≈ b_toml.structure.supercell.x_frac
    end

    @testset "AtomsBase path" begin
        sys = _periodic_pair()
        b = SCEBasis(sys;
            interaction = (body1 = (lmax = Dict(:Fe => 2),),
                           body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
        @test b isa SCEBasis
        @test b.structure.kd_name == ["Fe"]
        @test b.structure.supercell.num_atoms == 2
    end

    @testset "atom_name sublabels become distinct kinds" begin
        sys = _periodic_pair(atom_names = (:Fe_4a, :Fe_8e))
        b = SCEBasis(sys;
            interaction = (body1 = (lmax = Dict(:Fe_4a => 2, :Fe_8e => 2),),
                           body2 = (lsum = 2,
                                    cutoff = Dict((:Fe_4a, :Fe_4a) => -1.0,
                                                  (:Fe_4a, :Fe_8e) => -1.0,
                                                  (:Fe_8e, :Fe_8e) => -1.0))),
            verbosity = false,
        )
        @test Set(b.structure.kd_name) == Set(["Fe_4a", "Fe_8e"])
    end

    @testset "element-only key fans out to sublabels" begin
        sys = _periodic_pair(atom_names = (:Fe_4a, :Fe_8e))
        # element-only :Fe key — fans out to both Fe_4a and Fe_8e
        b = SCEBasis(sys;
            interaction = (body1 = (lmax = Dict(:Fe => 2),),
                           body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
        @test Set(b.structure.kd_name) == Set(["Fe_4a", "Fe_8e"])
    end

    @testset "odd lmax / lsum emits a warning" begin
        sys = _periodic_pair()
        @test_logs (:warn,) match_mode = :any SCEBasis(sys;
            interaction = (body1 = (lmax = Dict(:Fe => 3),),
                           body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
        @test_logs (:warn,) match_mode = :any SCEBasis(sys;
            interaction = (body1 = (lmax = Dict(:Fe => 2),),
                           body2 = (lsum = 3, cutoff = Dict((:Fe, :Fe) => -1.0))),
            verbosity = false,
        )
    end
end
