#!/usr/bin/env julia
"""
Convert VASP output to extended XYZ (extxyz) format.

Usage (minimal):
    julia vasp2extxyz.jl --vasprun vasprun.xml

Usage (with magnetic data):
    julia vasp2extxyz.jl --vasprun vasprun.xml --oszicar OSZICAR --output out.xyz
"""

include(joinpath(@__DIR__, "../ExtXYZ.jl"))
include(joinpath(@__DIR__, "VaspParser.jl"))

using .ExtXYZWriter
using .VaspParser
using ArgParse
using Printf

# ── comment string ─────────────────────────────────────────────────────────────

function build_vasp_comment(d::VaspRunData)::String
    parts = String[]
    !isempty(d.version)        && push!(parts, "VASP$(d.version)")
    d.encut        !== nothing && push!(parts, "ENCUT=$(Int(round(d.encut)))")
    d.kpoints_mesh !== nothing && push!(parts, "KPOINTS=$(d.kpoints_mesh)")
    d.iconst       !== nothing && push!(parts, "ICONST=$(d.iconst)")
    d.lambda       !== nothing && push!(parts, @sprintf("LAMBDA=%.4g", d.lambda))
    return join(parts, ", ")
end

# Expand per-species RWIGS to a 1×N per-atom matrix using atomtype indices.
# atomtype_per_atom[i] is the 1-based species-type index for atom i, matching
# the order of RWIGS values in the INCAR.  This correctly handles cases where
# the same element appears as multiple species with distinct RWIGS values.
function _rwigs_per_atom(rwigs::Vector{Float64}, atomtype_per_atom::Vector{Int})::Matrix{Float64}
    out = Matrix{Float64}(undef, 1, length(atomtype_per_atom))
    for (i, t) in enumerate(atomtype_per_atom)
        out[1, i] = rwigs[t]
    end
    return out
end

# ── CLI ────────────────────────────────────────────────────────────────────────

function parse_args_vasp2extxyz()
    s = ArgParseSettings(
        description = """
Convert VASP output to extxyz format.

When only --vasprun is given, a standard extxyz with structure, forces,
stress, and energies is produced.  Adding --oszicar also writes per-atom
magnetic moments (MW_int and M_int) and the constraint field (lambda*MW_perp).
""",
    )
    @add_arg_table! s begin
        "--vasprun", "-v"
            help     = "Path to vasprun.xml (required)"
            arg_type = String
            required = true
        "--oszicar", "-s"
            help     = "Path to OSZICAR — enables magnetic moment and constraint field columns"
            arg_type = String
            default  = nothing
        "--output", "-o"
            help     = "Output filename (default: stdout)"
            arg_type = String
            default  = nothing
    end
    return ArgParse.parse_args(ARGS, s)
end

function main()
    args = parse_args_vasp2extxyz()

    # ── parse vasprun.xml ──────────────────────────────────────────────────
    @info "Parsing $(args["vasprun"])…"
    vd = parse_vasprun(args["vasprun"])

    # Fractional → Cartesian positions (lattice columns = lattice vectors)
    positions_cart = vd.lattice * vd.positions_frac  # 3×N

    # Additional per-atom properties written as extra columns in extxyz.
    # Each entry is "name" => ncols×N matrix (ncols=1 for scalars, 3 for vectors).
    extra_per_atom = Pair{String, Matrix{Float64}}[]

    # ── RWIGS as per-atom scalar ───────────────────────────────────────────
    if vd.rwigs !== nothing
        push!(extra_per_atom, "rwigs" => _rwigs_per_atom(vd.rwigs, vd.atomtype_per_atom))
    end

    # ── parse OSZICAR (optional) ───────────────────────────────────────────
    if args["oszicar"] !== nothing
        @info "Parsing $(args["oszicar"]) for magnetic data…"
        md = parse_oszicar_magdata(args["oszicar"], vd.m_constr, vd.num_atoms)
        push!(extra_per_atom, "MAGMOM_smoothed" => md.magmom_smoothed)
        push!(extra_per_atom, "magmom_raw"      => md.magmom_raw)
        push!(extra_per_atom, "constr_field"    => md.constr_field)
    end

    # ── assemble frame ─────────────────────────────────────────────────────
    frame = AtomFrame(
        num_atoms     = vd.num_atoms,
        lattice       = vd.lattice,
        pbc           = vd.pbc,
        species       = vd.species,
        positions     = positions_cart,
        forces        = vd.forces,
        energy_free   = vd.energy_free,
        energy_zero   = vd.energy_zero,
        stress        = vd.stress,
        extra_per_atom = extra_per_atom,
        comment       = build_vasp_comment(vd),
    )

    # ── write output ───────────────────────────────────────────────────────
    if args["output"] !== nothing
        outfile = endswith(args["output"], ".extxyz") ? args["output"] : args["output"] * ".extxyz"
        write_extxyz(outfile, frame)
        @info "Written to $(outfile)"
    else
        write_extxyz(stdout, frame)
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
