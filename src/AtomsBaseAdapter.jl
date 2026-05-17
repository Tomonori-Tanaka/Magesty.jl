"""
    AtomsBaseAdapter

Internal adapter converting an `AtomsBase.AbstractSystem` plus an
`interaction` NamedTuple (or raw keyword inputs) into the typed value
objects (`SystemSpec`, `InteractionSpec`, `SymmetryOptions`) consumed
downstream. All AtomsBase / Unitful usage is confined to this file.
"""
module AtomsBaseAdapter

using AtomsBase
using OffsetArrays
using Unitful

using ..InputSpecs:
    SystemSpec,
    InteractionSpec,
    SymmetryOptions,
    expand_pair_table,
    expand_species_table

export system_to_specs, kwargs_to_specs

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
    for i = 1:n
        kind_syms[i] = _kind_name(sys, i)
        elem_syms[i] = Symbol(AtomsBase.atomic_symbol(AtomsBase.species(sys, i)))
    end
    unique_kinds = unique(kind_syms)
    kd_name = String.(unique_kinds)
    kd_int_list = [findfirst(==(k), unique_kinds)::Int for k in kind_syms]
    kind_to_element = Dict{Symbol, Symbol}()
    for i = 1:n
        kind_to_element[kind_syms[i]] = elem_syms[i]
    end
    return kd_name, kd_int_list, kind_to_element
end

# --- structure section --------------------------------------------------

"""
    _lattice_and_positions(sys) -> (lattice_mat, positions_frac_mat)

Lattice as a 3x3 `Matrix{Float64}` (columns = lattice vectors, angstrom)
and fractional atomic positions as a 3xN `Matrix{Float64}`. AtomsBase
stores Cartesian positions, so they are converted via `inv(lattice)`.
"""
function _lattice_and_positions(sys::AtomsBase.AbstractSystem)
    cv = AtomsBase.cell_vectors(sys)
    lattice_cols = [collect(Float64, _to_angstrom(c) for c in v) for v in cv]
    lattice_mat = reduce(hcat, lattice_cols)
    n = length(sys)
    positions_mat = zeros(Float64, 3, n)
    for i = 1:n
        cart = Float64[_to_angstrom(c) for c in AtomsBase.position(sys, i)]
        positions_mat[:, i] = lattice_mat \ cart
    end
    return lattice_mat, positions_mat
end

# --- interaction NamedTuple -> string-keyed dict via element fan-out ----

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

# Convert the body1.lmax NamedTuple/Dict to a fully-resolved
# Dict{String, Int} keyed by every kd_name (with element fan-out).
function _body1_lmax_dict(
    lmax_in,
    kd_name::AbstractVector{<:AbstractString},
    kind_to_element::AbstractDict{Symbol, Symbol},
)
    out = Dict{String, Int}()
    for kind in kd_name
        ksym = Symbol(kind)
        v = _lookup_species(lmax_in, ksym, kind_to_element[ksym], "lmax")
        out[kind] = _warn_if_odd(v, "lmax for species :$ksym")
    end
    return out
end

# Convert the body-n cutoff NamedTuple/Dict to a fully-resolved
# Dict{String, Float64} keyed by every unordered pair (with element fan-out).
function _bodyn_cutoff_dict(
    cutoff_in,
    kd_name::AbstractVector{<:AbstractString},
    kind_to_element::AbstractDict{Symbol, Symbol},
)
    out = Dict{String, Float64}()
    for i in eachindex(kd_name), j in i:length(kd_name)
        ki, kj = kd_name[i], kd_name[j]
        ksi, ksj = Symbol(ki), Symbol(kj)
        val = _lookup_pair(cutoff_in, ksi, ksj,
                           kind_to_element[ksi], kind_to_element[ksj])
        out["$ki-$kj"] = _to_angstrom(val)
    end
    return out
end

"""
    _build_interaction_spec(interaction, kd_name, kind_to_element) -> InteractionSpec

Build an `InteractionSpec` from the nested `interaction` NamedTuple.
Performs element-only key fan-out, cutoff unit normalization, and odd
`lmax` / `lsum` warnings, then constructs the typed spec.
"""
function _build_interaction_spec(
    interaction::NamedTuple,
    kd_name::AbstractVector{<:AbstractString},
    kind_to_element::AbstractDict{Symbol, Symbol},
)
    nbody = length(interaction)
    nbody >= 1 || throw(ArgumentError("interaction must have at least body1"))
    haskey(interaction, :body1) ||
        throw(ArgumentError("interaction must contain `body1`"))

    body1_lmax_dict = _body1_lmax_dict(interaction.body1.lmax, kd_name, kind_to_element)
    body1_lmax = expand_species_table(kd_name, body1_lmax_dict;
                                      context = "interaction.body1.lmax")

    nkd = length(kd_name)
    if nbody == 1
        bodyn_lsum = OffsetArray(Int[], 2:1)
        bodyn_cutoff = OffsetArray(zeros(Float64, 0, nkd, nkd), 2:1, 1:nkd, 1:nkd)
    else
        bodyn_lsum = OffsetArray(fill(0, nbody - 1), 2:nbody)
        bodyn_cutoff = OffsetArray(zeros(Float64, nbody - 1, nkd, nkd),
                                   2:nbody, 1:nkd, 1:nkd)
        for n = 2:nbody
            bkey = Symbol("body$n")
            haskey(interaction, bkey) || throw(ArgumentError(
                "interaction must contain `$bkey` for nbody = $nbody"))
            bn = getproperty(interaction, bkey)
            bodyn_lsum[n] = _warn_if_odd(bn.lsum, "lsum for body$n")
            cutoff_dict = _bodyn_cutoff_dict(bn.cutoff, kd_name, kind_to_element)
            cutoff_mat = expand_pair_table(kd_name, cutoff_dict;
                                           context = "interaction.body$n.cutoff")
            for i = 1:nkd, j = 1:nkd
                bodyn_cutoff[n, i, j] = cutoff_mat[i, j]
            end
        end
    end

    return InteractionSpec(
        nbody = nbody,
        body1_lmax = body1_lmax,
        bodyn_lsum = bodyn_lsum,
        bodyn_cutoff = bodyn_cutoff,
        kd_name = kd_name,
    )
end

# --- top-level ----------------------------------------------------------

"""
    system_to_specs(sys, interaction; name, tolerance_sym, isotropy)
        -> (SystemSpec, InteractionSpec, SymmetryOptions)

Build the three typed input value objects from an
`AtomsBase.AbstractSystem` and a nested `interaction` NamedTuple.

# Arguments
- `sys::AtomsBase.AbstractSystem`: The material system.
- `interaction::NamedTuple`: Nested `body1` / `body2` / ... interaction
  specification (see the spec's "Interaction parameters"). Element /
  kind fan-out is applied.
- `name::AbstractString = "system"`: System name recorded in `SystemSpec`.
- `tolerance_sym::Real = 1e-3`: Symmetry-detection tolerance.
- `isotropy::Bool = false`: Restrict to isotropic (Lf = 0) terms.
"""
function system_to_specs(
    sys::AtomsBase.AbstractSystem,
    interaction::NamedTuple;
    name::AbstractString = "system",
    tolerance_sym::Real = 1e-3,
    isotropy::Bool = false,
)
    kd_name, kd_int_list, kind_to_element = _kind_tables(sys)
    lattice_mat, positions_mat = _lattice_and_positions(sys)
    periodic = collect(Bool, AtomsBase.periodicity(sys))

    system = SystemSpec(
        name = String(name),
        num_atoms = length(sys),
        kd_name = kd_name,
        kd_int_list = kd_int_list,
        lattice_vectors = lattice_mat,
        x_fractional = positions_mat,
        is_periodic = periodic,
    )
    inter = _build_interaction_spec(interaction, kd_name, kind_to_element)
    options = SymmetryOptions(tolerance_sym = tolerance_sym, isotropy = isotropy)
    return (system, inter, options)
end

"""
    kwargs_to_specs(; lattice, kd, kd_list, positions, periodicity,
                       interaction, name, tolerance_sym, isotropy)
        -> (SystemSpec, InteractionSpec, SymmetryOptions)

Build the three typed input value objects from raw Julia keyword
arguments (no AtomsBase / Unitful required). `lattice` columns are the
lattice vectors (bare reals taken as angstrom, `Unitful.Length` also
accepted); `positions` are fractional coordinates; `kd` are the kind
names.

Unlike the AtomsBase path, the keyword path has no element/sublabel
distinction: every `interaction` key must match a `kd` entry exactly —
there is no element-only fan-out, since there is no element information
beyond `kd` itself.
"""
function kwargs_to_specs(;
    lattice::AbstractMatrix,
    kd::AbstractVector{Symbol},
    kd_list::AbstractVector{<:Integer},
    positions::AbstractVector,
    periodicity = (true, true, true),
    interaction::NamedTuple,
    name::AbstractString = "system",
    tolerance_sym::Real = 1e-3,
    isotropy::Bool = false,
)
    size(lattice) == (3, 3) ||
        throw(ArgumentError("lattice must be a 3x3 matrix, got $(size(lattice))"))
    kd_name = String.(kd)
    kind_to_element = Dict(Symbol(k) => Symbol(k) for k in kd_name)
    lattice_mat = Matrix{Float64}(undef, 3, 3)
    for i = 1:3, j = 1:3
        lattice_mat[i, j] = _to_angstrom(lattice[i, j])
    end
    n = length(kd_list)
    positions_mat = zeros(Float64, 3, n)
    for (i, p) in enumerate(positions)
        positions_mat[:, i] = collect(Float64, p)
    end
    periodic = collect(Bool, periodicity)

    system = SystemSpec(
        name = String(name),
        num_atoms = n,
        kd_name = kd_name,
        kd_int_list = collect(Int, kd_list),
        lattice_vectors = lattice_mat,
        x_fractional = positions_mat,
        is_periodic = periodic,
    )
    inter = _build_interaction_spec(interaction, kd_name, kind_to_element)
    options = SymmetryOptions(tolerance_sym = tolerance_sym, isotropy = isotropy)
    return (system, inter, options)
end

end # module AtomsBaseAdapter
