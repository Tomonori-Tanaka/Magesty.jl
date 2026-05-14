"""
    AtomsBaseAdapter

Internal adapter converting an `AtomsBase.AbstractSystem` plus an
`interaction` NamedTuple into the TOML-shaped input dictionary that
`ConfigParser.Config4System` already parses. This keeps all AtomsBase /
Unitful usage confined to one file and reuses the existing
config-parsing and validation logic unchanged.
"""
module AtomsBaseAdapter

using AtomsBase
using Unitful

export system_to_input_dict, kwargs_to_input_dict

# --- unit normalization -------------------------------------------------

# Distance-valued user input: a bare real is assumed to already be in
# angstrom, a Unitful length is converted. The -1.0 no-cutoff sentinel
# passes through the `Real` method unchanged.
_to_angstrom(x::Unitful.Length)::Float64 = ustrip(u"Å", x)
_to_angstrom(x::Real)::Float64 = Float64(x)

# --- per-atom kind name -------------------------------------------------

"""
    _kind_name(sys, i) -> Symbol

Kind name of atom `i`: the per-atom `:atom_name` data entry if present
(free-form, e.g. `:Fe_4a`), otherwise the element symbol. AtomsBase's
`ChemicalSpecies.atom_name` is deliberately not used — it is capped at
4 characters and ignored by `==`.
"""
function _kind_name(sys::AtomsBase.AbstractSystem, i::Integer)::Symbol
    atom = sys[i]
    if hasproperty(atom, :data) && haskey(atom.data, :atom_name)
        return Symbol(atom.data[:atom_name])
    end
    return Symbol(AtomsBase.atomic_symbol(AtomsBase.species(sys, i)))
end

"""
    _kind_tables(sys) -> (kd_name, kd_int_list, kind_to_element)

Build the unique kind-name list (`Vector{String}`), the per-atom kind
index (`Vector{Int}`), and a `kind Symbol -> element Symbol` map. Kind
order follows first appearance.
"""
function _kind_tables(sys::AtomsBase.AbstractSystem)
    n = length(sys)
    kind_syms = Vector{Symbol}(undef, n)
    elem_syms = Vector{Symbol}(undef, n)
    for i in 1:n
        kind_syms[i] = _kind_name(sys, i)
        elem_syms[i] = Symbol(AtomsBase.atomic_symbol(AtomsBase.species(sys, i)))
    end
    unique_kinds = unique(kind_syms)
    kd_name = String.(unique_kinds)
    kd_int_list = [findfirst(==(k), unique_kinds)::Int for k in kind_syms]
    kind_to_element = Dict{Symbol, Symbol}()
    for i in 1:n
        kind_to_element[kind_syms[i]] = elem_syms[i]
    end
    return kd_name, kd_int_list, kind_to_element
end

# --- structure section --------------------------------------------------

"""
    _lattice_and_positions(sys) -> (lattice_cols, positions_frac)

Lattice vectors as a `Vector` of three length-3 `Float64` vectors
(angstrom, each entry one lattice vector) and fractional atomic
positions as a `Vector` of length-3 `Float64` vectors. AtomsBase stores
Cartesian positions, so they are converted via `inv(lattice)`.
"""
function _lattice_and_positions(sys::AtomsBase.AbstractSystem)
    cv = AtomsBase.cell_vectors(sys)
    lattice_cols = [collect(Float64, _to_angstrom(c) for c in v) for v in cv]
    lattice_mat = reduce(hcat, lattice_cols)   # 3x3, columns = lattice vectors
    n = length(sys)
    positions_frac = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        cart = Float64[_to_angstrom(c) for c in AtomsBase.position(sys, i)]
        positions_frac[i] = lattice_mat \ cart
    end
    return lattice_cols, positions_frac
end

# --- interaction section ------------------------------------------------

# Look up an lmax-style per-species value: the kind name first, then the
# element symbol it belongs to (element-only key fan-out).
function _lookup_species(d::AbstractDict, kind::Symbol, element::Symbol,
                         what::AbstractString)
    haskey(d, kind) && return d[kind]
    haskey(d, element) && return d[element]
    throw(ArgumentError("$what not specified for kind :$kind (element :$element)"))
end

# Look up a cutoff value for a kind pair: the kind pair (either order)
# first, then the element pair (either order).
function _lookup_pair(d::AbstractDict, ki::Symbol, kj::Symbol,
                      ei::Symbol, ej::Symbol)
    for key in ((ki, kj), (kj, ki), (ei, ej), (ej, ei))
        haskey(d, key) && return d[key]
    end
    throw(ArgumentError(
        "cutoff not specified for kind pair (:$ki, :$kj) / " *
        "element pair (:$ei, :$ej)"))
end

function _warn_if_odd(value::Integer, label::AbstractString)
    if isodd(value)
        @warn "$label = $value is odd; time-reversal symmetry drops " *
              "odd-l terms, so this is equivalent to $label = $(value - 1)."
    end
    return value
end

"""
    _interaction_section(interaction, kd_name, kind_to_element) -> Dict

Convert the nested `interaction` NamedTuple into the TOML-shaped
`[interaction]` dictionary. Performs element-only key fan-out, cutoff
unit normalization, and odd `lmax` / `lsum` warnings.
"""
function _interaction_section(
    interaction::NamedTuple,
    kd_name::AbstractVector{<:AbstractString},
    kind_to_element::AbstractDict{Symbol, Symbol},
)
    nbody = length(interaction)
    section = Dict{String, Any}("nbody" => nbody)

    haskey(interaction, :body1) ||
        throw(ArgumentError("interaction must contain `body1`"))
    lmax_in = interaction.body1.lmax
    lmax_out = Dict{String, Int}()
    for kind in kd_name
        ksym = Symbol(kind)
        v = _lookup_species(lmax_in, ksym, kind_to_element[ksym], "lmax")
        lmax_out[kind] = _warn_if_odd(v, "lmax for species :$ksym")
    end
    section["body1"] = Dict{String, Any}("lmax" => lmax_out)

    for n in 2:nbody
        bkey = Symbol("body$n")
        haskey(interaction, bkey) || throw(ArgumentError(
            "interaction must contain `$bkey` for nbody = $nbody"))
        bn = getproperty(interaction, bkey)
        lsum = _warn_if_odd(bn.lsum, "lsum for body$n")
        cutoff_in = bn.cutoff
        cutoff_out = Dict{String, Float64}()
        for i in eachindex(kd_name), j in i:length(kd_name)
            ki, kj = kd_name[i], kd_name[j]
            ksi, ksj = Symbol(ki), Symbol(kj)
            val = _lookup_pair(cutoff_in, ksi, ksj,
                               kind_to_element[ksi], kind_to_element[ksj])
            cutoff_out["$ki-$kj"] = _to_angstrom(val)
        end
        section["body$n"] = Dict{String, Any}("lsum" => lsum, "cutoff" => cutoff_out)
    end

    return section
end

# --- top-level ----------------------------------------------------------

"""
    system_to_input_dict(sys, interaction; name, tolerance_sym, isotropy) -> Dict

Build the TOML-shaped input dictionary (`[general]` / `[symmetry]` /
`[interaction]` / `[structure]` sections) that `Config4System` parses,
from an `AtomsBase.AbstractSystem` and an `interaction` NamedTuple.

# Arguments
- `sys::AtomsBase.AbstractSystem`: The material system.
- `interaction::NamedTuple`: Nested `body1` / `body2` / ... interaction
  specification (see the spec's "Interaction parameters").
- `name::AbstractString = "system"`: System name recorded in `[general]`.
- `tolerance_sym::Real = 1e-5`: Symmetry-detection tolerance.
- `isotropy::Bool = false`: Restrict to isotropic (Lf = 0) terms.

# Returns
- `Dict{String, Any}`: The input dictionary, ready for `Config4System`.
"""
function system_to_input_dict(
    sys::AtomsBase.AbstractSystem,
    interaction::NamedTuple;
    name::AbstractString = "system",
    tolerance_sym::Real = 1e-5,
    isotropy::Bool = false,
)
    kd_name, kd_int_list, kind_to_element = _kind_tables(sys)
    lattice_cols, positions_frac = _lattice_and_positions(sys)
    periodic = collect(Bool, AtomsBase.periodicity(sys))

    return Dict{String, Any}(
        "general" => Dict{String, Any}(
            "name" => String(name),
            "nat" => length(sys),
            "kd" => kd_name,
            "periodicity" => periodic,
        ),
        "symmetry" => Dict{String, Any}(
            "tolerance" => Float64(tolerance_sym),
            "isotropy" => isotropy,
        ),
        "interaction" => _interaction_section(interaction, kd_name, kind_to_element),
        "structure" => Dict{String, Any}(
            "lattice" => lattice_cols,
            "kd_list" => kd_int_list,
            "position" => positions_frac,
        ),
    )
end

"""
    kwargs_to_input_dict(; lattice, kd, kd_list, positions, periodicity,
                         interaction, name, tolerance_sym, isotropy) -> Dict

Build the TOML-shaped input dictionary from raw Julia keyword arguments
(no AtomsBase / Unitful required). `lattice` columns are the lattice
vectors (bare reals taken as angstrom, `Unitful.Length` also accepted);
`positions` are fractional coordinates; `kd` are the kind names.

Unlike the AtomsBase path, the keyword path has no element/sublabel
distinction: every `interaction` key must match a `kd` entry exactly —
there is no element-only fan-out, since there is no element information
beyond `kd` itself.
"""
function kwargs_to_input_dict(;
    lattice::AbstractMatrix,
    kd::AbstractVector{Symbol},
    kd_list::AbstractVector{<:Integer},
    positions::AbstractVector,
    periodicity = (true, true, true),
    interaction::NamedTuple,
    name::AbstractString = "system",
    tolerance_sym::Real = 1e-5,
    isotropy::Bool = false,
)
    size(lattice) == (3, 3) ||
        throw(ArgumentError("lattice must be a 3x3 matrix, got $(size(lattice))"))
    kd_name = String.(kd)
    # keyword path: kind == element (no element-only fan-out)
    kind_to_element = Dict(Symbol(k) => Symbol(k) for k in kd_name)
    lattice_cols = [collect(Float64, _to_angstrom(lattice[i, j]) for i in 1:3)
                    for j in 1:3]
    positions_frac = [collect(Float64, p) for p in positions]
    periodic = collect(Bool, periodicity)

    return Dict{String, Any}(
        "general" => Dict{String, Any}(
            "name" => String(name),
            "nat" => length(kd_list),
            "kd" => kd_name,
            "periodicity" => periodic,
        ),
        "symmetry" => Dict{String, Any}(
            "tolerance" => Float64(tolerance_sym),
            "isotropy" => isotropy,
        ),
        "interaction" => _interaction_section(interaction, kd_name, kind_to_element),
        "structure" => Dict{String, Any}(
            "lattice" => lattice_cols,
            "kd_list" => collect(Int, kd_list),
            "position" => positions_frac,
        ),
    )
end

end # module AtomsBaseAdapter
