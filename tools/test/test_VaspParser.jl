"""
Integration tests for VaspParser, using fixture files under tools/test/fixtures/.

Run from repo root:
    julia tools/test/test_VaspParser.jl
"""

include(joinpath(@__DIR__, "../vasp/VaspParser.jl"))
using .VaspParser
using Test

const FIXTURE_DIR = joinpath(@__DIR__, "fixtures")
const FERH_DIR    = joinpath(FIXTURE_DIR, "FeRh")
const IRMN3_DIR   = joinpath(FIXTURE_DIR, "IrMn3")

# 1 kBar = 1/1602.1766208 eV/Å³  (same constant as VaspParser)
const KBAR_TO_EV_A3 = 1.0 / 1602.1766208

# ── parse_vasprun: FeRh ───────────────────────────────────────────────────────

@testset "parse_vasprun: FeRh" begin
    vd = parse_vasprun(joinpath(FERH_DIR, "vasprun.xml"))

    @test vd.num_atoms == 16
    @test vd.version   == "6.6.0"

    # Species: 8 Fe then 8 Rh
    @test vd.species[1]  == "Fe"
    @test vd.species[8]  == "Fe"
    @test vd.species[9]  == "Rh"
    @test vd.species[16] == "Rh"

    # Atomtype indices: Fe→1, Rh→2
    @test all(vd.atomtype_per_atom[1:8]  .== 1)
    @test all(vd.atomtype_per_atom[9:16] .== 2)

    # INCAR scalars
    @test vd.encut  ≈ 300.0
    @test vd.lambda ≈ 2.0
    @test vd.iconst == 1

    # RWIGS (one value per species)
    @test vd.rwigs !== nothing
    @test length(vd.rwigs) == 2
    @test vd.rwigs[1] ≈ 1.24470
    @test vd.rwigs[2] ≈ 1.34030

    # k-mesh
    @test vd.kpoints_mesh == "3 3 3"

    # Lattice: cubic 5.97 Å (columns = lattice vectors)
    @test vd.lattice[1, 1] ≈ 5.97
    @test vd.lattice[2, 2] ≈ 5.97
    @test vd.lattice[3, 3] ≈ 5.97
    @test vd.lattice[1, 2] ≈ 0.0 atol=1e-10

    # Periodic boundary conditions
    @test vd.pbc == [true, true, true]

    # Fractional positions: atom 1 at origin, atom 2 at (0, 0, 0.5)
    @test all(vd.positions_frac[:, 1] .≈ 0.0)
    @test vd.positions_frac[3, 2] ≈ 0.5

    # Energies (last calculation block)
    @test vd.energy_free ≈ -124.80455729 atol=1e-6
    @test vd.energy_zero ≈ -124.79753671 atol=1e-6

    # Stress (last calculation block); VASP output in kBar → eV/Å³
    @test vd.stress[1, 1] ≈ 33.25821644 * KBAR_TO_EV_A3 atol=1e-8
    @test vd.stress[3, 3] ≈ 33.42401004 * KBAR_TO_EV_A3 atol=1e-8

    # M_CONSTR: 3×16, all non-zero (Fe z=3, Rh z=1)
    @test vd.m_constr !== nothing
    @test size(vd.m_constr) == (3, 16)
    @test vd.m_constr[3, 1] ≈ 3.0    # Fe: M_CONSTR z-component
    @test vd.m_constr[3, 9] ≈ 1.0    # Rh: M_CONSTR z-component
end

# ── parse_vasprun: IrMn3 ─────────────────────────────────────────────────────

@testset "parse_vasprun: IrMn3" begin
    vd = parse_vasprun(joinpath(IRMN3_DIR, "vasprun.xml"))

    @test vd.num_atoms == 4

    # Species: Mn×3, Ir×1
    @test vd.species[1] == "Mn"
    @test vd.species[3] == "Mn"
    @test vd.species[4] == "Ir"

    # Atomtype: Mn→1, Ir→2
    @test all(vd.atomtype_per_atom[1:3] .== 1)
    @test vd.atomtype_per_atom[4] == 2

    # INCAR
    @test vd.encut  ≈ 350.0
    @test vd.lambda ≈ 1.0
    @test vd.iconst == 4

    # RWIGS
    @test vd.rwigs !== nothing
    @test vd.rwigs[1] ≈ 1.302
    @test vd.rwigs[2] ≈ 1.455

    # k-mesh
    @test vd.kpoints_mesh == "5 5 5"

    # Energies
    @test vd.energy_free ≈ -36.71506424 atol=1e-6
    @test vd.energy_zero ≈ -36.71579408 atol=1e-6

    # M_CONSTR: last atom (Ir) is all-zero → unconstrained
    @test vd.m_constr !== nothing
    @test size(vd.m_constr) == (3, 4)
    @test all(abs.(vd.m_constr[:, 4]) .< 1e-10)   # Ir: unconstrained
    @test vd.m_constr[1, 1] ≈ 3.0                  # Mn1 x
end

# ── parse_oszicar_magdata: FeRh ───────────────────────────────────────────────

@testset "parse_oszicar_magdata: FeRh (all atoms constrained)" begin
    vd = parse_vasprun(joinpath(FERH_DIR, "vasprun.xml"))
    md = parse_oszicar_magdata(joinpath(FERH_DIR, "OSZICAR"), vd.m_constr, vd.num_atoms)

    @test size(md.magmom_smoothed) == (3, 16)
    @test size(md.magmom_raw)      == (3, 16)
    @test size(md.constr_field)    == (3, 16)

    # Last step: Fe (atoms 1–8) MW_int z ≈ 2.189
    for i in 1:8
        @test md.magmom_smoothed[3, i] ≈ 2.189 atol=1e-3
    end

    # Last step: Rh (atoms 9–16) MW_int z ≈ 0.663
    for i in 9:16
        @test md.magmom_smoothed[3, i] ≈ 0.663 atol=1e-3
    end

    # All atoms are constrained → constr_field should be non-zero (or at least
    # populated); actual values are ~1e-18 (negligible but non-zero in OSZICAR)
    # We check the structure, not the tiny numerical values
    @test !all(iszero, md.constr_field)   # at least some entries read

    # x and y components for Fe should be ≈ 0 (only z constrained)
    @test md.magmom_smoothed[1, 1] ≈ 0.0 atol=1e-3
    @test md.magmom_smoothed[2, 1] ≈ 0.0 atol=1e-3
end

# ── parse_oszicar_magdata: IrMn3 (last atom unconstrained) ───────────────────

@testset "parse_oszicar_magdata: IrMn3 (Ir atom unconstrained)" begin
    vd = parse_vasprun(joinpath(IRMN3_DIR, "vasprun.xml"))
    md = parse_oszicar_magdata(joinpath(IRMN3_DIR, "OSZICAR"), vd.m_constr, vd.num_atoms)

    @test size(md.magmom_smoothed) == (3, 4)

    # Last step MW_int values (from OSZICAR tail)
    @test md.magmom_smoothed[1, 1] ≈  1.484 atol=1e-3
    @test md.magmom_smoothed[2, 1] ≈  0.742 atol=1e-3
    @test md.magmom_smoothed[3, 1] ≈ -0.743 atol=1e-3

    @test md.magmom_smoothed[1, 2] ≈ -0.742 atol=1e-3
    @test md.magmom_smoothed[2, 2] ≈ -1.485 atol=1e-3
    @test md.magmom_smoothed[3, 2] ≈ -0.743 atol=1e-3

    @test md.magmom_smoothed[1, 3] ≈ -0.742 atol=1e-3
    @test md.magmom_smoothed[2, 3] ≈  0.742 atol=1e-3
    @test md.magmom_smoothed[3, 3] ≈  1.484 atol=1e-3

    # lambda*MW_perp: Ir (atom 4, M_CONSTR=0) → all zeros
    @test all(md.constr_field[:, 4] .== 0.0)

    # Mn atoms (1–3) should have nonzero constraint field
    @test md.constr_field[1, 1] ≈ -0.37197e-4 atol=1e-8
    @test md.constr_field[2, 1] ≈ -1.1692e-4  atol=1e-8
    @test md.constr_field[3, 1] ≈ -1.9125e-4  atol=1e-8

    @test md.constr_field[1, 2] ≈  0.23545e-3 atol=1e-8
    @test md.constr_field[1, 3] ≈  0.15883e-4 atol=1e-8
end

println("All VaspParser tests passed.")
