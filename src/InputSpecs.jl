module InputSpecs

using OffsetArrays

export SystemSpec, InteractionSpec, SymmetryOptions, parse_toml_inputs

const DEFAULT_TOLERANCE_SYM = 1e-3
const DEFAULT_ISOTROPY = false
const DEFAULT_PERIODICITY = [true, true, true]

"""
    SystemSpec

Material-side input: name, atom count, species, lattice, fractional
positions, periodicity. Holds no interaction or symmetry parameters.

# Fields
- `name::String`: System identifier.
- `num_atoms::Int`: Number of atoms in the cell.
- `kd_name::Vector{String}`: Unique species labels in user-defined order.
- `kd_int_list::Vector{Int}`: Per-atom species index (1-based, into `kd_name`).
- `lattice_vectors::Matrix{Float64}`: 3 x 3 matrix; columns are a, b, c (Angstrom).
- `x_fractional::Matrix{Float64}`: 3 x num_atoms fractional coordinates.
- `is_periodic::Vector{Bool}`: Length 3, per-direction periodicity.
"""
struct SystemSpec
    name::String
    num_atoms::Int
    kd_name::Vector{String}
    kd_int_list::Vector{Int}
    lattice_vectors::Matrix{Float64}
    x_fractional::Matrix{Float64}
    is_periodic::Vector{Bool}

    function SystemSpec(;
        name::AbstractString,
        num_atoms::Integer,
        kd_name::AbstractVector{<:AbstractString},
        kd_int_list::AbstractVector{<:Integer},
        lattice_vectors::AbstractMatrix{<:Real},
        x_fractional::AbstractMatrix{<:Real},
        is_periodic::AbstractVector{Bool} = DEFAULT_PERIODICITY,
    )
        isempty(name) && throw(ArgumentError("Structure name cannot be empty"))
        num_atoms > 0 ||
            throw(ArgumentError("Number of atoms must be positive, got $num_atoms"))
        isempty(kd_name) &&
            throw(ArgumentError("Chemical species list cannot be empty"))
        allunique(kd_name) ||
            throw(ArgumentError("Chemical species list must be unique, got $(kd_name)"))
        size(lattice_vectors) == (3, 3) ||
            throw(ArgumentError("Lattice vectors must be a 3x3 matrix, got $(size(lattice_vectors))"))
        size(x_fractional) == (3, num_atoms) ||
            throw(ArgumentError("x_fractional must be 3x$num_atoms, got $(size(x_fractional))"))
        length(kd_int_list) == num_atoms ||
            throw(ArgumentError("kd_int_list length ($(length(kd_int_list))) must match num_atoms ($num_atoms)"))
        nkd = length(kd_name)
        for (i, k) in pairs(kd_int_list)
            (1 <= k <= nkd) || throw(
                ArgumentError(
                    "kd_int_list[$i] = $k out of range 1:$nkd",
                ),
            )
        end
        length(is_periodic) == 3 ||
            throw(ArgumentError("Periodicity must be specified for all three directions, got length $(length(is_periodic))"))
        return new(
            String(name),
            Int(num_atoms),
            Vector{String}(kd_name),
            Vector{Int}(kd_int_list),
            Matrix{Float64}(lattice_vectors),
            Matrix{Float64}(x_fractional),
            Vector{Bool}(is_periodic),
        )
    end
end

"""
    InteractionSpec

Interaction-side input: many-body order and per-body / per-species
angular momentum and cutoff parameters.

# Fields
- `nbody::Int`: Maximum interaction body order (>= 1).
- `body1_lmax::Vector{Int}`: Body-1 lmax per species, length = `length(kd_name)`.
- `bodyn_lsum::OffsetArray{Int, 1}`: lsum per body, indexed `2:nbody`.
  Empty (`2:1`) when `nbody == 1`.
- `bodyn_cutoff::OffsetArray{Float64, 3}`: Pair cutoff per (body, kd_i, kd_j),
  indexed `(2:nbody, 1:nkd, 1:nkd)`. Symmetric in the last two indices.
  `-1.0` means "include all pairs". Empty body axis (`2:1`) when `nbody == 1`.

# Construction
`kd_name::Vector{String}` is consumed for cross-validation against the
species labels but not stored.
"""
struct InteractionSpec
    nbody::Int
    body1_lmax::Vector{Int}
    bodyn_lsum::OffsetArray{Int, 1}
    bodyn_cutoff::OffsetArray{Float64, 3}

    function InteractionSpec(;
        nbody::Integer,
        body1_lmax::AbstractVector{<:Integer},
        bodyn_lsum::AbstractArray{<:Integer, 1},
        bodyn_cutoff::AbstractArray{<:Real, 3},
        kd_name::AbstractVector{<:AbstractString},
    )
        nbody > 0 ||
            throw(ArgumentError("nbody must be positive, got $nbody"))
        nkd = length(kd_name)
        length(body1_lmax) == nkd || throw(
            ArgumentError(
                "body1_lmax length ($(length(body1_lmax))) must match number of species ($nkd)",
            ),
        )
        if nbody == 1
            length(bodyn_lsum) == 0 || throw(
                ArgumentError(
                    "bodyn_lsum must be empty when nbody == 1, got length $(length(bodyn_lsum))",
                ),
            )
            size(bodyn_cutoff, 1) == 0 || throw(
                ArgumentError(
                    "bodyn_cutoff body axis must be empty when nbody == 1",
                ),
            )
        else
            axes(bodyn_lsum) == (2:nbody,) || throw(
                ArgumentError(
                    "bodyn_lsum axes must be (2:$nbody,), got $(axes(bodyn_lsum))",
                ),
            )
            axes(bodyn_cutoff) == (2:nbody, 1:nkd, 1:nkd) || throw(
                ArgumentError(
                    "bodyn_cutoff axes must be (2:$nbody, 1:$nkd, 1:$nkd), got $(axes(bodyn_cutoff))",
                ),
            )
            for n in 2:nbody
                for i in 1:nkd, j in 1:nkd
                    bodyn_cutoff[n, i, j] == bodyn_cutoff[n, j, i] || throw(
                        ArgumentError(
                            "bodyn_cutoff[$n, $i, $j] = $(bodyn_cutoff[n, i, j]) not symmetric with [$n, $j, $i] = $(bodyn_cutoff[n, j, i])",
                        ),
                    )
                end
            end
        end
        return new(
            Int(nbody),
            Vector{Int}(body1_lmax),
            OffsetArray(Int.(parent(bodyn_lsum)), axes(bodyn_lsum)),
            OffsetArray(Float64.(parent(bodyn_cutoff)), axes(bodyn_cutoff)),
        )
    end
end

"""
    SymmetryOptions(; tolerance_sym = 1e-3, isotropy = false)

Symmetry-detection options.

# Fields
- `tolerance_sym::Float64`: Tolerance for spglib symmetry detection (> 0).
- `isotropy::Bool`: If `true`, keep only isotropic basis terms (Lf = 0).
"""
struct SymmetryOptions
    tolerance_sym::Float64
    isotropy::Bool

    function SymmetryOptions(;
        tolerance_sym::Real = DEFAULT_TOLERANCE_SYM,
        isotropy::Bool = DEFAULT_ISOTROPY,
    )
        tolerance_sym > 0 || throw(
            ArgumentError("tolerance_sym must be positive, got $tolerance_sym"),
        )
        return new(Float64(tolerance_sym), isotropy)
    end
end


# --------------------------------------------------------------------------
# Wildcard-aware parsing helpers
# --------------------------------------------------------------------------

# Internal representation of one pair-table entry after key parsing.
struct _PairEntry
    key::String
    tier::Int                            # 0 = concrete-concrete, 1 = one wildcard, 2 = both wildcard
    a::Union{Int, Nothing}               # species index or nothing for '*'
    b::Union{Int, Nothing}
    value::Float64
end

# Internal representation of one species-table entry.
struct _SpeciesEntry
    key::String
    tier::Int                            # 0 = concrete, 1 = wildcard
    idx::Union{Int, Nothing}
    value::Int
end

function _parse_pair_key(key::AbstractString, kd_name::AbstractVector{<:AbstractString})
    parts = split(key, "-")
    length(parts) == 2 || throw(
        ArgumentError("malformed pair key \"$key\" (expected \"<a>-<b>\")"),
    )
    a, b = String.(parts)
    a_star = (a == "*")
    b_star = (b == "*")
    a_idx = a_star ? nothing : findfirst(==(a), kd_name)
    b_idx = b_star ? nothing : findfirst(==(b), kd_name)
    if !a_star && a_idx === nothing
        throw(ArgumentError("unknown species \"$a\" in pair key \"$key\""))
    end
    if !b_star && b_idx === nothing
        throw(ArgumentError("unknown species \"$b\" in pair key \"$key\""))
    end
    tier = Int(a_star) + Int(b_star)
    return (a_idx = a_idx, b_idx = b_idx, tier = tier)
end

# Canonical form for duplicate detection. Two distinct input keys collapse to
# the same canonical iff they specify the same (unordered) pair pattern.
function _canonical_pair(info, kd_name::AbstractVector{<:AbstractString})
    if info.tier == 2
        return "*-*"
    elseif info.tier == 1
        concrete = info.a_idx === nothing ? kd_name[info.b_idx] : kd_name[info.a_idx]
        return string(concrete, "-*")
    else
        s1 = kd_name[info.a_idx]
        s2 = kd_name[info.b_idx]
        return s1 <= s2 ? string(s1, "-", s2) : string(s2, "-", s1)
    end
end

function _matches_pair(p::_PairEntry, i::Int, j::Int)
    a_ok = (p.a === nothing) || (p.a == i)
    b_ok = (p.b === nothing) || (p.b == j)
    if a_ok && b_ok
        return true
    end
    a_ok2 = (p.a === nothing) || (p.a == j)
    b_ok2 = (p.b === nothing) || (p.b == i)
    return a_ok2 && b_ok2
end

"""
    expand_pair_table(kd_name, entries; context = "")

Expand a pair table that may contain wildcard keys (`"*-*"`, `"X-*"`,
`"*-Y"`) into a fully concrete symmetric `Matrix{Float64}` of size
`(nkd, nkd)`, where `nkd = length(kd_name)`. Resolution is
specificity-based: concrete-concrete > one wildcard > all wildcard.

Throws `ArgumentError` for duplicate equivalent keys, unknown species,
ambiguous same-tier coverage with different values, or missing coverage
of any pair.
"""
function expand_pair_table(
    kd_name::AbstractVector{<:AbstractString},
    entries::AbstractDict;
    context::AbstractString = "",
)::Matrix{Float64}
    nkd = length(kd_name)
    parsed = _PairEntry[]
    canonical_seen = Dict{String, String}()
    for (raw_key, raw_value) in entries
        key = String(raw_key)
        info = _parse_pair_key(key, kd_name)
        canon = _canonical_pair(info, kd_name)
        if haskey(canonical_seen, canon)
            prev = canonical_seen[canon]
            throw(ArgumentError(
                "duplicate pair specification: \"$prev\" and \"$key\" both designate the same pair$(_ctx(context))",
            ))
        end
        canonical_seen[canon] = key
        push!(parsed, _PairEntry(key, info.tier, info.a_idx, info.b_idx, Float64(raw_value)))
    end

    out = zeros(Float64, nkd, nkd)
    for i in 1:nkd, j in i:nkd
        matches = [p for p in parsed if _matches_pair(p, i, j)]
        if isempty(matches)
            throw(ArgumentError(
                "no cutoff specified for pair \"$(kd_name[i])-$(kd_name[j])\"$(_ctx(context))",
            ))
        end
        best_tier = minimum(p.tier for p in matches)
        winners = [p for p in matches if p.tier == best_tier]
        if length(winners) > 1
            distinct_vals = unique(p.value for p in winners)
            if length(distinct_vals) > 1
                keys_str = join((p.key for p in winners), ", ")
                throw(ArgumentError(
                    "ambiguous wildcard coverage for pair \"$(kd_name[i])-$(kd_name[j])\": keys [$keys_str] tie at tier $best_tier with conflicting values$(_ctx(context))",
                ))
            end
        end
        v = winners[1].value
        out[i, j] = v
        out[j, i] = v
    end
    return out
end

"""
    expand_species_table(kd_name, entries; context = "")

Expand a per-species table that may contain a wildcard key (`"*"`) into
a concrete `Vector{Int}` of length `nkd`. Concrete species keys override
the wildcard for that species.

Returns `zeros(Int, nkd)` if `entries` is empty (preserves the
"body1 section omitted" default of the legacy parser).
"""
function expand_species_table(
    kd_name::AbstractVector{<:AbstractString},
    entries::AbstractDict;
    context::AbstractString = "",
)::Vector{Int}
    nkd = length(kd_name)
    if isempty(entries)
        return zeros(Int, nkd)
    end
    parsed = _SpeciesEntry[]
    for (raw_key, raw_value) in entries
        key = String(raw_key)
        if key == "*"
            push!(parsed, _SpeciesEntry(key, 1, nothing, Int(raw_value)))
        else
            idx = findfirst(==(key), kd_name)
            idx === nothing && throw(
                ArgumentError("unknown species \"$key\" in species table$(_ctx(context))"),
            )
            push!(parsed, _SpeciesEntry(key, 0, idx, Int(raw_value)))
        end
    end
    # tier-1 (wildcard) keys must be unique by construction (Dict has one "*").
    out = fill(typemin(Int), nkd)
    has_wildcard = any(e -> e.tier == 1, parsed)
    wildcard_val = has_wildcard ? first(e for e in parsed if e.tier == 1).value : 0
    for i in 1:nkd
        concrete_hit = filter(e -> e.tier == 0 && e.idx == i, parsed)
        if !isempty(concrete_hit)
            out[i] = concrete_hit[1].value
        elseif has_wildcard
            out[i] = wildcard_val
        else
            throw(ArgumentError(
                "no value specified for species \"$(kd_name[i])\"$(_ctx(context))",
            ))
        end
    end
    return out
end

_ctx(s::AbstractString) = isempty(s) ? "" : " ($s)"

# --------------------------------------------------------------------------
# parse_toml_inputs: single dict/TOML entry point
# --------------------------------------------------------------------------

const _REQUIRED_SECTIONS = ("general", "symmetry", "interaction", "structure")

"""
    parse_toml_inputs(dict::AbstractDict) -> (SystemSpec, InteractionSpec, SymmetryOptions)

Parse a TOML-shape dictionary into the three typed input value objects.
Wildcard species notation (`"*-*"`, `"X-*"`, `"*-Y"` for pair tables and
`"*"` for species tables) is expanded here using specificity-based
resolution; the returned `InteractionSpec` is fully concrete (no
wildcards survive).

Expected schema (matches the legacy `Config4System` schema; no breaking
changes):

```toml
[general]
name        = "..."
nat         = N
kd          = ["A", "B", ...]
periodicity = [true, true, true]    # optional, default [true, true, true]

[symmetry]
tolerance = 1e-3                    # optional, default 1e-3
isotropy  = false                   # optional, default false

[interaction]
nbody = M
[interaction.body1.lmax]
"A" = ...
"*" = ...                           # newly accepted

[interaction.body2]
lsum = ...
[interaction.body2.cutoff]
"A-A" = ...
"A-*" = ...                         # newly accepted
"*-*" = ...                         # newly accepted

[structure]
lattice  = [[...], [...], [...]]
kd_list  = [1, ..., M]
position = [[x, y, z], ...]
```
"""
function parse_toml_inputs(dict::AbstractDict)
    for section in _REQUIRED_SECTIONS
        haskey(dict, section) || throw(
            ArgumentError("Required section \"$section\" is missing in the input dictionary."),
        )
    end

    general = dict["general"]
    sym     = dict["symmetry"]
    inter   = dict["interaction"]
    str_    = dict["structure"]

    kd_name = Vector{String}(general["kd"])
    nat = Int(general["nat"])

    lattice_vectors = hcat([Float64.(v) for v in str_["lattice"]]...)
    positions_raw = [Float64.(vec) for vec in str_["position"]]
    x_fractional = _parse_positions(positions_raw, nat)
    kd_int_list = Vector{Int}(str_["kd_list"])
    is_periodic = Vector{Bool}(get(general, "periodicity", DEFAULT_PERIODICITY))

    system = SystemSpec(
        name = String(general["name"]),
        num_atoms = nat,
        kd_name = kd_name,
        kd_int_list = kd_int_list,
        lattice_vectors = lattice_vectors,
        x_fractional = x_fractional,
        is_periodic = is_periodic,
    )

    nbody = Int(inter["nbody"])
    body1_lmax = if haskey(inter, "body1") && haskey(inter["body1"], "lmax")
        expand_species_table(kd_name, inter["body1"]["lmax"]; context = "interaction.body1.lmax")
    else
        zeros(Int, length(kd_name))
    end

    nkd = length(kd_name)
    if nbody == 1
        bodyn_lsum = OffsetArray(Int[], 2:1)
        bodyn_cutoff = OffsetArray(zeros(Float64, 0, nkd, nkd), 2:1, 1:nkd, 1:nkd)
    else
        bodyn_lsum = OffsetArray(fill(typemin(Int), nbody - 1), 2:nbody)
        bodyn_cutoff = OffsetArray(zeros(Float64, nbody - 1, nkd, nkd), 2:nbody, 1:nkd, 1:nkd)
        for n in 2:nbody
            bodyn_key = "body$n"
            haskey(inter, bodyn_key) || throw(
                ArgumentError("interaction.$bodyn_key section is missing for nbody = $nbody"),
            )
            bodyn = inter[bodyn_key]
            haskey(bodyn, "lsum") || throw(
                ArgumentError("interaction.$bodyn_key.lsum is missing"),
            )
            bodyn_lsum[n] = Int(bodyn["lsum"])
            haskey(bodyn, "cutoff") || throw(
                ArgumentError("interaction.$bodyn_key.cutoff is missing"),
            )
            cutoff_mat = expand_pair_table(
                kd_name, bodyn["cutoff"];
                context = "interaction.$bodyn_key.cutoff",
            )
            for i in 1:nkd, j in 1:nkd
                bodyn_cutoff[n, i, j] = cutoff_mat[i, j]
            end
        end
    end

    interaction = InteractionSpec(
        nbody = nbody,
        body1_lmax = body1_lmax,
        bodyn_lsum = bodyn_lsum,
        bodyn_cutoff = bodyn_cutoff,
        kd_name = kd_name,
    )

    options = SymmetryOptions(
        tolerance_sym = get(sym, "tolerance", DEFAULT_TOLERANCE_SYM),
        isotropy = get(sym, "isotropy", DEFAULT_ISOTROPY),
    )

    return (system, interaction, options)
end

function _parse_positions(
    position_list::AbstractVector{<:AbstractVector{<:Real}},
    num_atoms::Integer,
)::Matrix{Float64}
    if isempty(position_list)
        throw(ArgumentError("Position list cannot be empty"))
    end
    if length(position_list) != num_atoms
        throw(ArgumentError(
            "Number of positions ($(length(position_list))) must match number of atoms ($num_atoms)",
        ))
    end
    out = fill(0.0, 3, num_atoms)
    for (i, vec) in enumerate(position_list)
        length(vec) == 3 || throw(ArgumentError(
            "Position vector for atom $i must have 3 components, got $(length(vec)).",
        ))
        out[:, i] = vec
    end
    return out
end

end # module
