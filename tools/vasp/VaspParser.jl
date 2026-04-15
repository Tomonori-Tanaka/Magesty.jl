"""
    VaspParser

Parse VASP output files (vasprun.xml and OSZICAR) into plain Julia structs.
"""
module VaspParser

using EzXML
using Printf

export VaspRunData, OszicarMagData, parse_vasprun, parse_oszicar_magdata

# 1 kBar = 0.1 GPa; 1 eV/Å³ ≈ 1602.18 GPa  →  1 kBar ≈ 6.2415e-4 eV/Å³
const KBAR_TO_EV_A3 = 1.0 / 1602.1766208

"""
Data extracted from a single (last) ionic step in vasprun.xml.

Fields
- `version`       : VASP version string (empty string if absent)
- `encut`         : ENCUT (eV), or `nothing`
- `kpoints_mesh`  : k-mesh string e.g. "4 4 4", or `nothing`
- `iconst`        : I_CONSTRAINED_M (or ICONST), or `nothing`
- `lambda`        : LAMBDA (magnetic penalty), or `nothing`
- `rwigs`         : RWIGS vector (one value per species), or `nothing`
- `m_constr`      : 3×N constraint direction matrix, or `nothing`
- `lattice`       : 3×3, columns = lattice vectors (Å)
- `species`            : element symbol per atom, length N
- `atomtype_per_atom`  : 1-based species-type index per atom (from atominfo); used to map RWIGS correctly
- `positions_frac`: 3×N fractional coordinates
- `forces`        : 3×N forces (eV/Å)
- `stress`        : 3×3 stress tensor (eV/Å³)
- `energy_free`   : free energy F (eV)
- `energy_zero`   : energy sigma→0 (eV)
- `num_atoms`     : N
- `pbc`           : always [true, true, true] for VASP
"""
struct VaspRunData
    version::String
    encut::Union{Float64, Nothing}
    kpoints_mesh::Union{String, Nothing}
    iconst::Union{Int, Nothing}
    lambda::Union{Float64, Nothing}
    rwigs::Union{Vector{Float64}, Nothing}
    m_constr::Union{Matrix{Float64}, Nothing}  # 3×N
    lattice::Matrix{Float64}                   # 3×3, columns = lattice vectors
    species::Vector{String}
    atomtype_per_atom::Vector{Int}             # 1-based species-type index per atom
    positions_frac::Matrix{Float64}            # 3×N
    forces::Matrix{Float64}                    # 3×N (eV/Å)
    stress::Matrix{Float64}                    # 3×3 (eV/Å³)
    energy_free::Float64
    energy_zero::Float64
    num_atoms::Int
    pbc::Vector{Bool}
end

"""
Magnetic data extracted from the last complete SCF step in OSZICAR.

Fields
- `magmom_smoothed`: 3×N Wannier-interpolated magnetic moments (MW_int in VASP)
- `magmom_raw`     : 3×N directly integrated magnetic moments (M_int in VASP)
- `constr_field`   : 3×N constraint (penalty) field λ·MW_perp; zero for unconstrained atoms
"""
struct OszicarMagData
    magmom_smoothed::Matrix{Float64}
    magmom_raw::Matrix{Float64}
    constr_field::Matrix{Float64}
end

# ── vasprun.xml ────────────────────────────────────────────────────────────────

"""
    parse_vasprun(filepath) -> VaspRunData

Parse the last ionic step from a vasprun.xml file.
"""
function parse_vasprun(filepath::AbstractString)::VaspRunData
    isfile(filepath) || throw(ArgumentError("File not found: $filepath"))
    doc = readxml(filepath)

    version      = _read_version(doc)
    incar_node   = findfirst("//incar", doc)
    encut        = _incar_scalar(incar_node, "ENCUT",  Float64)
    lambda       = _incar_scalar(incar_node, "LAMBDA", Float64)
    # Try both tag names used across VASP versions for the constraint type
    iconst       = something(_incar_scalar(incar_node, "I_CONSTRAINED_M", Int),
                             _incar_scalar(incar_node, "ICONST", Int))
    rwigs        = _incar_vector(incar_node, "RWIGS")
    m_constr_flat = _incar_vector(incar_node, "M_CONSTR")

    kpoints_mesh      = _read_kpoints_mesh(doc)
    species, atomtype_per_atom = _read_atominfo(doc)
    num_atoms         = length(species)

    # Use the last <calculation> block
    calc_nodes = findall("//calculation", doc)
    isempty(calc_nodes) && error("No <calculation> element found in $filepath")
    calc = calc_nodes[end]

    lattice, positions_frac = _read_structure(calc)
    forces  = _read_varray(calc, "forces")
    stress  = _read_stress(calc)
    energy_free, energy_zero = _read_energies(calc)

    m_constr = nothing
    if m_constr_flat !== nothing && length(m_constr_flat) == 3 * num_atoms
        m_constr = reshape(m_constr_flat, 3, num_atoms)
    end

    return VaspRunData(
        version, encut, kpoints_mesh, iconst, lambda, rwigs, m_constr,
        lattice, species, atomtype_per_atom, positions_frac, forces, stress,
        energy_free, energy_zero, num_atoms, [true, true, true],
    )
end

function _read_version(doc)::String
    node = findfirst("//generator/i[@name='version']", doc)
    node === nothing && return ""
    return strip(nodecontent(node))
end

function _incar_scalar(incar_node, name::AbstractString, ::Type{T}) where {T}
    incar_node === nothing && return nothing
    node = findfirst("i[@name='$(name)']", incar_node)
    node === nothing && return nothing
    try
        return parse(T, strip(nodecontent(node)))
    catch
        return nothing
    end
end

function _incar_vector(incar_node, name::AbstractString)::Union{Vector{Float64}, Nothing}
    incar_node === nothing && return nothing
    # Try single <v name="..."> first (most common for vector INCAR tags)
    node = findfirst("v[@name='$(name)']", incar_node)
    if node !== nothing
        vals = parse.(Float64, split(strip(nodecontent(node))))
        return isempty(vals) ? nothing : vals
    end
    # Fallback: <varray name="..."><v>...</v></varray>
    arr_node = findfirst("varray[@name='$(name)']", incar_node)
    arr_node === nothing && return nothing
    vals = Float64[]
    for v in findall("v", arr_node)
        append!(vals, parse.(Float64, split(strip(nodecontent(v)))))
    end
    return isempty(vals) ? nothing : vals
end

function _read_kpoints_mesh(doc)::Union{String, Nothing}
    gen_node = findfirst("//kpoints/generation", doc)
    gen_node === nothing && return nothing
    div_node = findfirst("v[@name='divisions']", gen_node)
    div_node === nothing && return nothing
    vals = parse.(Int, split(strip(nodecontent(div_node))))
    return join(string.(vals), " ")
end

# Returns (species, atomtype_per_atom).
# Each <rc> in the atoms array has two <c> children: element name and 1-based atomtype index.
# Using the atomtype index (not the element name) to map RWIGS correctly handles
# the case where the same element appears as multiple species with different RWIGS values.
function _read_atominfo(doc)::Tuple{Vector{String}, Vector{Int}}
    rc_nodes = findall("//atominfo/array[@name='atoms']/set/rc", doc)
    isempty(rc_nodes) && error("Cannot find atom species in vasprun.xml")
    species   = String[]
    atomtypes = Int[]
    for rc in rc_nodes
        cs = findall("c", rc)
        length(cs) < 2 && error("Unexpected <rc> format in atominfo")
        push!(species,   strip(nodecontent(cs[1])))
        push!(atomtypes, parse(Int, strip(nodecontent(cs[2]))))
    end
    return species, atomtypes
end

function _read_structure(calc_node)
    struct_node = findfirst("structure", calc_node)
    struct_node === nothing && error("No <structure> in calculation node")

    basis_nodes = findall("crystal/varray[@name='basis']/v", struct_node)
    isempty(basis_nodes) && error("No basis vectors found in <structure>")
    # Each <v> is one lattice vector; hcat → columns = lattice vectors
    lattice = hcat([parse.(Float64, split(strip(nodecontent(v)))) for v in basis_nodes]...)

    pos_nodes = findall("varray[@name='positions']/v", struct_node)
    positions_frac = hcat([parse.(Float64, split(strip(nodecontent(v)))) for v in pos_nodes]...)

    return lattice, positions_frac
end

function _read_varray(calc_node, name::AbstractString)::Matrix{Float64}
    arr_node = findfirst("varray[@name='$(name)']", calc_node)
    arr_node === nothing && error("No <varray name=\"$(name)\"> in calculation")
    rows = [parse.(Float64, split(strip(nodecontent(v)))) for v in findall("v", arr_node)]
    return hcat(rows...)  # 3×N
end

function _read_stress(calc_node)::Matrix{Float64}
    arr_node = findfirst("varray[@name='stress']", calc_node)
    arr_node === nothing && return zeros(3, 3)
    rows = [parse.(Float64, split(strip(nodecontent(v)))) for v in findall("v", arr_node)]
    stress_kbar = hcat(rows...)'   # 3×3, row i = i-th stress row
    return Matrix(stress_kbar) .* KBAR_TO_EV_A3
end

function _read_energies(calc_node)
    e_node = findfirst("energy", calc_node)
    e_node === nothing && error("No <energy> in last calculation")
    e_free = parse(Float64, nodecontent(findfirst("i[@name='e_fr_energy']", e_node)))
    e_zero = parse(Float64, nodecontent(findfirst("i[@name='e_0_energy']",  e_node)))
    return e_free, e_zero
end

# ── OSZICAR ────────────────────────────────────────────────────────────────────

"""
    parse_oszicar_magdata(filepath, m_constr, num_atoms) -> OszicarMagData

Parse magnetic moment and constraint field from the *last* complete SCF step
in an OSZICAR (or OUTCAR) file.

The magnetic moment section is identified by a header line containing
"ion", "MW_int", and "M_int", followed by rows of 7 columns:
    ion_idx  MWx  MWy  MWz  Mx  My  Mz

The constraint field section starts with a line containing "lambda*MW_perp",
followed by rows of 4 columns: ion_idx  cfx  cfy  cfz
VASP omits atoms for which M_CONSTR = 0 0 0; those are set to zero here.
"""
function parse_oszicar_magdata(
    filepath::AbstractString,
    m_constr::Union{Matrix{Float64}, Nothing},
    num_atoms::Int,
)::OszicarMagData
    isfile(filepath) || throw(ArgumentError("File not found: $filepath"))

    # Atoms with M_CONSTR == 0 will not appear in the lambda*MW_perp section
    is_constrained = fill(true, num_atoms)
    if m_constr !== nothing
        for i in 1:num_atoms
            if all(abs(m_constr[k, i]) < 1e-10 for k in 1:3)
                is_constrained[i] = false
            end
        end
    end

    # Working buffers; committed to result arrays whenever a section ends
    magmom_smoothed = zeros(3, num_atoms)
    magmom_raw      = zeros(3, num_atoms)
    constr_field    = zeros(3, num_atoms)
    tmp_smoothed = zeros(3, num_atoms)
    tmp_raw      = zeros(3, num_atoms)
    tmp_cf       = zeros(3, num_atoms)

    in_magmom = false
    in_constr = false

    for line in eachline(filepath)
        stripped = strip(line)

        # ── magmom section ──────────────────────────────────────────────────
        if !in_magmom && occursin("ion", stripped) &&
                occursin("MW_int", stripped) && occursin("M_int", stripped)
            in_magmom = true
            fill!(tmp_smoothed, 0.0)
            fill!(tmp_raw,      0.0)
            continue
        end

        if in_magmom
            parts = split(stripped)
            if length(parts) == 7
                try
                    idx = parse(Int, parts[1])
                    if 1 <= idx <= num_atoms
                        tmp_smoothed[:, idx] .= parse.(Float64, parts[2:4])
                        tmp_raw[:, idx]      .= parse.(Float64, parts[5:7])
                    end
                    continue
                catch; end
            end
            # Section ended — commit and fall through to check other headers
            in_magmom = false
            magmom_smoothed .= tmp_smoothed
            magmom_raw      .= tmp_raw
        end

        # ── constraint field section ────────────────────────────────────────
        if !in_constr && occursin("lambda*MW_perp", stripped)
            in_constr = true
            fill!(tmp_cf, 0.0)
            continue
        end

        if in_constr
            parts = split(stripped)
            if length(parts) == 4
                try
                    idx = parse(Int, parts[1])
                    if 1 <= idx <= num_atoms && is_constrained[idx]
                        tmp_cf[:, idx] .= parse.(Float64, parts[2:4])
                    end
                    continue
                catch; end
            end
            # Section ended — commit
            in_constr = false
            constr_field .= tmp_cf
        end
    end

    # Capture sections that extend to the end of file
    if in_magmom
        magmom_smoothed .= tmp_smoothed
        magmom_raw      .= tmp_raw
    end
    if in_constr
        constr_field .= tmp_cf
    end

    return OszicarMagData(magmom_smoothed, magmom_raw, constr_field)
end

end # module VaspParser
