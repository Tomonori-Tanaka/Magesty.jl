# Tool for mean-field analysis of SCE models.
# Reads jphi.xml / scecoeffs.xml, Fourier-transforms exchange interactions J(q),
# and finds the ordering wave vector with minimum MFA energy (= maximum eigenvalue
# of J(q)).
using ArgParse
using LinearAlgebra
using Printf

using Magesty
include("convert2tensor.jl")
using .ExchangeTensor

"""
    prim_lattice_in_sc_frac(symmetry) -> Matrix{Float64}

Return the 3×3 matrix whose columns are primitive cell lattice vectors expressed
in supercell fractional coordinates, constructed from the pure translation
symmetry operations.  Falls back to the identity if the system is already
primitive (no non-trivial translations found).
"""
function prim_lattice_in_sc_frac(symmetry)::Matrix{Float64}
    tol = 1e-5
    vecs = [Vector{Float64}(symmetry.symdata[i].translation_frac)
            for i in symmetry.symnum_translation
            if norm(symmetry.symdata[i].translation_frac) > tol]

    isempty(vecs) && return Matrix{Float64}(I, 3, 3)
    sort!(vecs, by = norm)

    # Pick the 3 shortest linearly independent translation vectors
    basis = Vector{Vector{Float64}}()
    for t in vecs
        length(basis) >= 3 && break
        M = isempty(basis) ? reshape(t, 3, 1) : hcat(basis..., t)
        if rank(M) > length(basis)
            push!(basis, t)
        end
    end

    length(basis) < 3 && return Matrix{Float64}(I, 3, 3)
    return hcat(basis...)
end

"""
    build_jq_matrix(q_sc, symmetry, structure, pre) -> Matrix{ComplexF64}

Construct the Fourier-transformed exchange matrix J(q) at wavevector `q_sc`
given in supercell fractional reciprocal coordinates.

J(q) is a (3 n_sub × 3 n_sub) block matrix.  Block (μ, ν) is:

    J_μν(q) = Σ_{atom_j ∈ sublattice ν} J(prim_μ, atom_j) exp(2πi q·R_{μj})

where R_{μj} is the minimum-image displacement (supercell fractional coords)
from the primitive-cell representative of sublattice μ to atom_j, and
J(prim_μ, atom_j) is the 3×3 exchange tensor in meV.
"""
function build_jq_matrix(
    q_sc::AbstractVector{<:Real},
    symmetry,
    structure,
    pre::ExchangeTensor.PrecomputedXMLData,
)::Matrix{ComplexF64}
    atoms_in_prim = symmetry.atoms_in_prim
    n_sub = length(atoms_in_prim)
    x_frac = structure.supercell.x_frac
    num_atoms = structure.supercell.num_atoms

    J_q = zeros(ComplexF64, 3n_sub, 3n_sub)

    for (μ, prim_μ) in enumerate(atoms_in_prim)
        rows = (3(μ-1)+1):(3μ)
        for atom_j in 1:num_atoms
            atom_j == prim_μ && continue  # skip on-site (R=0, same atom)

            # Sublattice index of atom_j (index into atoms_in_prim, 1-based)
            ν = symmetry.map_s2p[atom_j].atom
            cols = (3(ν-1)+1):(3ν)

            # Minimum-image displacement in supercell fractional coordinates
            R = x_frac[:, atom_j] .- x_frac[:, prim_μ]
            R .= R .- round.(R)

            phase = exp(2π * im * dot(q_sc, R))
            tensor = ExchangeTensor.convert2tensor_fast(symmetry, prim_μ, atom_j, pre)
            J_q[rows, cols] .+= tensor.jij_tensor .* phase
        end
    end

    return J_q
end

"""
    mfa_analysis(xml_path; nk, spin) -> (q_opt, λ_max, eigenvalues)

Scan `nk^3` q-points uniformly over the primitive cell first Brillouin zone
and return the ordering wave vector that minimizes the smallest eigenvalue of
J(q), corresponding to the minimum MFA energy.

Arguments:
- `xml_path`: path to jphi.xml or scecoeffs.xml
- `nk`      : number of q-points per reciprocal direction (default 20)
- `spin`    : spin magnitude S.
              If `nothing` (default): classical spin, Tc ∝ S²/3  with S=1.
              If specified: quantum spin, Tc ∝ S(S+1)/3.

Returns `(q_opt_prim, λ_min, eigenvalues_at_qopt)` where `q_opt_prim` is in
primitive cell fractional reciprocal coordinates.
"""
function mfa_analysis(xml_path::String; nk::Int = 20, spin::Union{Float64, Nothing} = nothing)
    println("Loading from: $xml_path")
    structure = Magesty.Structure(xml_path, verbosity = false)
    symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)
    pre = ExchangeTensor.precompute_xml(xml_path)

    atoms_in_prim = symmetry.atoms_in_prim
    n_sub = length(atoms_in_prim)
    println("Primitive cell atoms (sublattices): $atoms_in_prim  [n_sub=$n_sub]")

    A_prim = prim_lattice_in_sc_frac(symmetry)
    println("Primitive lattice vectors in supercell frac coords (columns of A_prim):")
    for col in eachcol(A_prim)
        @printf("  [%+.4f  %+.4f  %+.4f]\n", col...)
    end
    @printf("Scanning %d^3 = %d q-points in the primitive BZ...\n", nk, nk^3)

    # Scan q_frac_prim ∈ [0, 1)^3, convert to supercell fractional for the phase
    # Convention in Magesty: H = +Σ J_ij S_i·S_j with J < 0 for ferromagnetic coupling.
    # The MFA ground state minimizes E(q) ∝ λ_min(J(q)), so we seek the q that gives
    # the most negative minimum eigenvalue of J(q).
    grid = range(0.0, 1.0, length = nk + 1)[1:nk]
    λ_min = +Inf
    q_prim_best = zeros(3)
    eigs_best = Float64[]

    for qx in grid, qy in grid, qz in grid
        q_prim = [qx, qy, qz]
        q_sc = A_prim * q_prim

        J_q = build_jq_matrix(q_sc, symmetry, structure, pre)
        J_q_herm = Hermitian((J_q + J_q') / 2)
        eigs = real(eigvals(J_q_herm))
        λ = minimum(eigs)

        if λ < λ_min
            λ_min = λ
            q_prim_best = q_prim
            eigs_best = sort(eigs)
        end
    end

    # Compute wavelength from q in Cartesian reciprocal space (1/Å)
    # A_prim_cart: columns = primitive lattice vectors in Å
    # B_prim = inv(A_prim_cart): rows = b_i reciprocal vectors (b_i · a_j = δ_ij, 1/Å)
    # q_cart = B_prim' * q_prim_best  (1/Å)
    # wavelength = 1 / |q_cart|  (Å)
    A_prim_cart = Matrix(structure.supercell.lattice_vectors) * A_prim
    B_prim = inv(A_prim_cart)
    q_cart = B_prim' * q_prim_best
    q_norm = norm(q_cart)

    println("\n=== Mean-Field Analysis Results ===")
    @printf("Optimal ordering vector q (primitive reciprocal frac): [%.4f, %.4f, %.4f]\n",
        q_prim_best...)
    @printf("Ordering vector |q| (Cartesian):  %.6f Å⁻¹\n", q_norm)
    if q_norm < 1e-10
        println("Wavelength: ∞  (ferromagnetic, q = 0)")
    else
        wavelength = 1.0 / q_norm
        a1_norm = norm(A_prim_cart[:, 1])
        @printf("Wavelength:   %.4f Å  =  %.4f a_prim\n", wavelength, wavelength / a1_norm)
    end
    @printf("Min eigenvalue λ_min of J(q):  %+.4f meV\n", λ_min)
    if isnothing(spin)
        # Classical spin: Tc = S²/3 * |λ_min|  with S=1
        Tc = 1.0 / 3 * (-λ_min) / 8.617333e-2
        @printf("Estimated MFA Tc (classical):  %.1f K\n", Tc)
    else
        # Quantum spin: Tc = S(S+1)/3 * |λ_min|
        Tc = spin * (spin + 1) / 3 * (-λ_min) / 8.617333e-2
        @printf("Estimated MFA Tc (quantum, S=%.1f):  %.1f K\n", spin, Tc)
    end
    println("All eigenvalues at q_opt (meV):")
    for (i, ev) in enumerate(eigs_best)
        @printf("  λ_%d = %+.4f meV\n", i, ev)
    end

    return q_prim_best, λ_min, eigs_best
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--xml", "-x"
            help = "Path to jphi.xml or scecoeffs.xml"
            arg_type = String
            required = true
        "--nk", "-n"
            help = "Number of q-points per reciprocal direction"
            arg_type = Int
            default = 20
        "--spin", "-s"
            help = "Spin magnitude S. If omitted: classical spin (Tc ∝ S²/3, S=1). If given: quantum spin (Tc ∝ S(S+1)/3)."
            arg_type = Float64
            default = nothing
    end
    args = parse_args(s)
    mfa_analysis(args["xml"], nk = args["nk"], spin = args["spin"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
