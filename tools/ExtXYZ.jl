"""
    ExtXYZ

Generic writer for the extended XYZ (extxyz) format.
Separating this module from VASP-specific code allows future support for other codes.

Reference: https://github.com/libAtoms/extxyz
"""
module ExtXYZWriter

using Printf

export AtomFrame, write_extxyz

"""
A single frame of atomic data for extxyz output.

Fields
- `lattice`       : 3×3 matrix whose *columns* are lattice vectors (Å)
- `pbc`           : 3-element periodic boundary flags
- `species`       : element symbols, length N
- `positions`     : 3×N Cartesian coordinates (Å)
- `forces`        : 3×N forces (eV/Å), or `nothing`
- `energy_free`   : free energy (eV), or `nothing`
- `energy_zero`   : sigma→0 energy (eV), or `nothing`
- `stress`        : 3×3 stress tensor (eV/Å³), or `nothing`
- `soc`           : spin-orbit coupling flag (true/false), or `nothing`
- `code`          : DFT code name (e.g. "VASP", "QE"), or `nothing`
- `version`       : DFT code version string, or `nothing`
- `extra_per_atom`: ordered list of additional per-atom real properties as (name => matrix) pairs.
                    Each matrix is `ncols×N`; ncols=1 for scalars, ncols=3 for vectors.
                    The Properties descriptor becomes `name:R:ncols` automatically.
- `comment`       : string placed in the `comment=` key of the extxyz header line
"""
struct AtomFrame
    num_atoms::Int
    lattice::Matrix{Float64}
    pbc::Vector{Bool}
    species::Vector{String}
    positions::Matrix{Float64}
    forces::Union{Matrix{Float64}, Nothing}
    energy_free::Union{Float64, Nothing}
    energy_zero::Union{Float64, Nothing}
    stress::Union{Matrix{Float64}, Nothing}
    soc::Union{Bool, Nothing}
    code::Union{String, Nothing}
    version::Union{String, Nothing}
    extra_per_atom::Vector{Pair{String, Matrix{Float64}}}
    comment::String
end

function AtomFrame(;
    num_atoms::Int,
    lattice::Matrix{Float64},
    pbc::Vector{Bool} = [true, true, true],
    species::Vector{String},
    positions::Matrix{Float64},
    forces::Union{Matrix{Float64}, Nothing} = nothing,
    energy_free::Union{Float64, Nothing} = nothing,
    energy_zero::Union{Float64, Nothing} = nothing,
    stress::Union{Matrix{Float64}, Nothing} = nothing,
    soc::Union{Bool, Nothing} = nothing,
    code::Union{String, Nothing} = nothing,
    version::Union{String, Nothing} = nothing,
    extra_per_atom::Vector{Pair{String, Matrix{Float64}}} = Pair{String, Matrix{Float64}}[],
    comment::String = "",
)
    return AtomFrame(
        num_atoms, lattice, pbc, species, positions,
        forces, energy_free, energy_zero, stress, soc,
        code, version, extra_per_atom, comment,
    )
end

"""
    write_extxyz(io, frame)
    write_extxyz(filename, frame)

Write one `AtomFrame` in extxyz format.
"""
function write_extxyz(io::IO, frame::AtomFrame)
    println(io, frame.num_atoms)
    println(io, _header_line(frame))
    for i in 1:frame.num_atoms
        cols = String[frame.species[i]]
        for k in 1:3
            push!(cols, @sprintf("% .10f", frame.positions[k, i]))
        end
        if frame.forces !== nothing
            for k in 1:3
                push!(cols, @sprintf("% .10f", frame.forces[k, i]))
            end
        end
        for (_, mat) in frame.extra_per_atom
            for k in axes(mat, 1)
                push!(cols, @sprintf("% .10f", mat[k, i]))
            end
        end
        println(io, join(cols, " "))
    end
end

function write_extxyz(filename::AbstractString, frame::AtomFrame)
    open(filename, "w") do io
        write_extxyz(io, frame)
    end
end

# ── internals ──────────────────────────────────────────────────────────────────

function _header_line(frame::AtomFrame)::String
    # Lattice="a1x a1y a1z  a2x a2y a2z  a3x a3y a3z"
    # lattice[:,i] = i-th lattice vector, so iterate column-major for row-major output
    lat = join([@sprintf("%.10f", frame.lattice[j, i]) for i in 1:3 for j in 1:3], " ")

    # Properties string: each descriptor is "name:type:ncols"
    # ncols is inferred from the first dimension of the per-atom matrix
    prop_descs = String["species:S:1", "pos:R:3"]
    frame.forces !== nothing && push!(prop_descs, "forces:R:3")
    for (name, mat) in frame.extra_per_atom
        push!(prop_descs, "$(name):R:$(size(mat, 1))")
    end
    props = join(prop_descs, ":")

    pbc_str = join([b ? "T" : "F" for b in frame.pbc], " ")

    parts = String[
        "Lattice=\"$(lat)\"",
        "Properties=$(props)",
        "pbc=\"$(pbc_str)\"",
    ]

    frame.code        !== nothing && push!(parts, "code=$(_quote_if_needed(frame.code))")
    frame.version     !== nothing && push!(parts, "version=$(_quote_if_needed(frame.version))")
    frame.energy_free !== nothing && push!(parts, @sprintf("energy_free=%.10f", frame.energy_free))
    frame.energy_zero !== nothing && push!(parts, @sprintf("energy_zero=%.10f", frame.energy_zero))
    frame.soc !== nothing         && push!(parts, "soc=$(frame.soc ? "T" : "F")")

    if frame.stress !== nothing
        S = frame.stress
        # Voigt notation (eV/Å³): xx yy zz yz xz xy
        vs = @sprintf("%.6e %.6e %.6e %.6e %.6e %.6e",
            S[1,1], S[2,2], S[3,3], S[2,3], S[1,3], S[1,2])
        push!(parts, "stress=\"$(vs)\"")
    end

    !isempty(frame.comment) && push!(parts, "comment=\"$(frame.comment)\"")

    return join(parts, " ")
end

# Quote string values that contain whitespace; leave bare strings unquoted.
function _quote_if_needed(s::AbstractString)::String
    return occursin(r"\s", s) ? "\"$(s)\"" : String(s)
end

end # module ExtXYZWriter
