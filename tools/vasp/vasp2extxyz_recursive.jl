#!/usr/bin/env julia
"""
Recursively convert VASP output to extxyz format.

Walks all subdirectories under --root (default: current directory), looking for
vasprun.xml and optionally OSZICAR in each directory.

Output modes (--mode):
  each      Write <dirname>.extxyz inside each processed directory.
  combined  Append all frames to a single file (--output, default: combined.extxyz).
  both      Do both of the above.

Skipping / warning rules:
  - Neither vasprun.xml nor OSZICAR found → silently skip the directory.
  - Only one of the two files found → warn and proceed with what is available.
  - Parse error in either file → warn and skip that directory.
"""

include(joinpath(@__DIR__, "../ExtXYZ.jl"))
include(joinpath(@__DIR__, "VaspParser.jl"))

using .ExtXYZWriter
using .VaspParser
using ArgParse
using Printf

# ── shared helpers (same as vasp2extxyz.jl) ───────────────────────────────────

function build_vasp_comment(d::VaspRunData)::String
    parts = String[]
    !isempty(d.version)        && push!(parts, "VASP$(d.version)")
    d.encut        !== nothing && push!(parts, "ENCUT=$(Int(round(d.encut)))")
    d.kpoints_mesh !== nothing && push!(parts, "KPOINTS=$(d.kpoints_mesh)")
    d.iconst       !== nothing && push!(parts, "ICONST=$(d.iconst)")
    d.lambda       !== nothing && push!(parts, @sprintf("LAMBDA=%.4g", d.lambda))
    return join(parts, ", ")
end

function _rwigs_per_atom(rwigs::Vector{Float64}, atomtype_per_atom::Vector{Int})::Matrix{Float64}
    out = Matrix{Float64}(undef, 1, length(atomtype_per_atom))
    for (i, t) in enumerate(atomtype_per_atom)
        out[1, i] = rwigs[t]
    end
    return out
end

# ── frame builder ─────────────────────────────────────────────────────────────

"""
    build_frame(dir, vasprun_name, oszicar_name) -> AtomFrame

Parse VASP output in `dir` and return an `AtomFrame`.  `oszicar_name` may be
`nothing` to skip magnetic data.  Throws on parse failure.
"""
function build_frame(
    dir::AbstractString,
    vasprun_name::AbstractString,
    oszicar_name::Union{AbstractString, Nothing},
)::AtomFrame
    vd = parse_vasprun(joinpath(dir, vasprun_name))

    positions_cart = vd.lattice * vd.positions_frac

    extra_per_atom = Pair{String, Matrix{Float64}}[]

    if vd.rwigs !== nothing
        push!(extra_per_atom, "rwigs" => _rwigs_per_atom(vd.rwigs, vd.atomtype_per_atom))
    end

    if oszicar_name !== nothing
        md = parse_oszicar_magdata(joinpath(dir, oszicar_name), vd.m_constr, vd.num_atoms)
        push!(extra_per_atom, "MAGMOM_smoothed" => md.magmom_smoothed)
        push!(extra_per_atom, "magmom_raw"      => md.magmom_raw)
        push!(extra_per_atom, "constr_field"    => md.constr_field)
    end

    return AtomFrame(
        num_atoms      = vd.num_atoms,
        lattice        = vd.lattice,
        pbc            = vd.pbc,
        species        = vd.species,
        positions      = positions_cart,
        forces         = vd.forces,
        energy_free    = vd.energy_free,
        energy_zero    = vd.energy_zero,
        stress         = vd.stress,
        extra_per_atom = extra_per_atom,
        comment        = build_vasp_comment(vd),
    )
end

# ── CLI ───────────────────────────────────────────────────────────────────────

function parse_args_recursive()
    s = ArgParseSettings(
        description = """
Recursively convert VASP output directories to extxyz.

Output modes:
  each      Write <dirname>.extxyz inside each processed directory.
  combined  Append all frames to one file (see --output).
  both      Do both of the above.
""",
    )
    @add_arg_table! s begin
        "--root", "-r"
            help     = "Root directory to search (default: current directory)"
            arg_type = String
            default  = "."
        "--vasprun", "-v"
            help     = "vasprun.xml filename to look for (default: vasprun.xml)"
            arg_type = String
            default  = "vasprun.xml"
        "--oszicar", "-s"
            help     = "OSZICAR filename to look for (default: OSZICAR)"
            arg_type = String
            default  = "OSZICAR"
        "--mode", "-m"
            help     = "Output mode: each | combined | both (default: each)"
            arg_type = String
            default  = "each"
        "--output", "-o"
            help     = "Output filename for combined mode (default: combined.extxyz)"
            arg_type = String
            default  = "combined.extxyz"
    end
    args = ArgParse.parse_args(ARGS, s)
    args["mode"] in ("each", "combined", "both") ||
        error("--mode must be one of: each, combined, both")
    return args
end

# ── main ──────────────────────────────────────────────────────────────────────

function main()
    args = parse_args_recursive()

    root         = args["root"]
    vasprun_name = args["vasprun"]
    oszicar_name = args["oszicar"]
    mode         = args["mode"]
    combined_out = endswith(args["output"], ".extxyz") ? args["output"] : args["output"] * ".extxyz"

    isdir(root) || error("Root directory not found: $root")

    # Collect all directories (root itself + all descendants), sorted for
    # reproducible ordering.
    all_dirs = String[]
    for (dirpath, _, _) in walkdir(root)
        push!(all_dirs, dirpath)
    end
    sort!(all_dirs)

    combined_io = (mode in ("combined", "both")) ? open(combined_out, "w") : nothing

    n_ok      = 0
    n_skipped = 0
    n_errors  = 0

    for dir in all_dirs
        has_vasprun = isfile(joinpath(dir, vasprun_name))
        has_oszicar = isfile(joinpath(dir, oszicar_name))

        # Neither file → silent skip
        if !has_vasprun && !has_oszicar
            n_skipped += 1
            continue
        end

        # Only one file present → warn, still attempt conversion
        if !has_vasprun
            println("WARNING: $(dir) — $(oszicar_name) found but $(vasprun_name) missing; skipping")
            n_errors += 1
            continue
        end
        if !has_oszicar
            println("NOTE: $(dir) — $(vasprun_name) found, $(oszicar_name) not found; magnetic data will be omitted")
        end

        # Attempt to parse and build the frame
        frame = try
            build_frame(dir, vasprun_name, has_oszicar ? oszicar_name : nothing)
        catch e
            println("ERROR: $(dir) — $(sprint(showerror, e))")
            n_errors += 1
            continue
        end

        # Write outputs
        if mode in ("each", "both")
            each_path = joinpath(dir, basename(dir) * ".extxyz")
            try
                write_extxyz(each_path, frame)
            catch e
                println("ERROR: $(dir) — could not write $(each_path): $(sprint(showerror, e))")
                n_errors += 1
                continue
            end
        end

        if combined_io !== nothing
            write_extxyz(combined_io, frame)
        end

        n_ok += 1
    end

    combined_io !== nothing && close(combined_io)

    @info "Done: $(n_ok) converted, $(n_errors) with warnings/errors, $(n_skipped) skipped (no VASP files)"
    if mode in ("combined", "both")
        @info "Combined output: $(combined_out)"
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
