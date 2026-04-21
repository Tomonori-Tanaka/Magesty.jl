#!/usr/bin/env julia
# Benchmark: MySphericalHarmonics vs SpheriCart in the Optimize.jl hot path.
#
# Compares three variants for design_matrix_energy_element and calc_∇ₑu:
#   baseline      – current code: Zₗₘ_unsafe / ∂ᵢZlm_unsafe per (l,m) per atom
#   sphericart_v1 – SpheriCart.compute / compute_with_gradients per atom (drop-in)
#   sphericart_v2 – all harmonics precomputed once per spin-config (restructured)
#
# Run via:
#   make bench-optimize-sphericart
#
# ENV overrides:
#   BENCH_SAMPLES=10 BENCH_EVALS=3 make bench-optimize-sphericart

using BenchmarkTools
using Printf
using StaticArrays
using LinearAlgebra
using TOML
using Magesty
using Magesty.Optimize: design_matrix_energy_element, calc_∇ₑu,
                        build_design_matrix_energy, build_design_matrix_torque
using Magesty.MySphericalHarmonics: Zₗₘ_unsafe, ∂ᵢZlm_unsafe
using SpheriCart

const BENCH_SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "10"))
const BENCH_EVALS   = parse(Int, get(ENV, "BENCH_EVALS",   "3"))
const BENCH_EXAMPLE = get(ENV, "BENCH_EXAMPLE", "fege_2x2x2")

# SpheriCart stores harmonics as Y(0,0), Y(1,-1), Y(1,0), Y(1,1), Y(2,-2), ...
# For m_idx in 1:(2l+1), the SpheriCart index is l^2 + m_idx (1-based).
@inline _sph_idx(l, m_idx) = l * l + m_idx
_sph_lmax(::SphericalHarmonics{L}) where L = L  # L is a type parameter, not a field

# ── sphericart_v1: one compute() call per atom, same call structure ────────────

function design_matrix_energy_element_v1(
    cbc,
    spin_directions::AbstractMatrix{<:Real},
    symmetry,
    sph::SphericalHarmonics,
)::Float64
    result = 0.0
    N = length(cbc.atoms)
    dims = [2*l + 1 for l in cbc.ls]
    site_indices = CartesianIndices(Tuple(dims))
    Mf_size = size(cbc.coeff_tensor, N + 1)
    searched_pairs = Set{Vector{Int}}()

    for itrans in symmetry.symnum_translation
        translated_atoms = [symmetry.map_sym[atom, itrans] for atom in cbc.atoms]
        atoms_sorted = sort(translated_atoms)
        atoms_sorted in searched_pairs && continue
        push!(searched_pairs, atoms_sorted)

        sh_values = Vector{Vector{Float64}}(undef, N)
        for (site_idx, atom) in enumerate(translated_atoms)
            l = cbc.ls[site_idx]
            Y = compute(sph, SVector{3,Float64}(spin_directions[:, atom]...))
            sh_values[site_idx] = [Y[_sph_idx(l, mi)] for mi in 1:(2l+1)]
        end

        tensor_result = 0.0
        for mf_idx in 1:Mf_size
            mf_contribution = 0.0
            for site_idx_tuple in site_indices
                product = 1.0
                for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
                    product *= sh_values[site_idx][m_idx]
                end
                mf_contribution += cbc.coeff_tensor[site_idx_tuple.I..., mf_idx] * product
            end
            tensor_result += cbc.coefficient[mf_idx] * mf_contribution
        end
        result += tensor_result * cbc.multiplicity
    end
    return result
end

function calc_∇ₑu_v1(
    cbc,
    atom::Integer,
    spin_directions::AbstractMatrix{<:Real},
    symmetry,
    sph::SphericalHarmonics,
)::Vector{Float64}
    result = MVector{3,Float64}(0.0, 0.0, 0.0)
    N = length(cbc.atoms)
    dims = [2*l + 1 for l in cbc.ls]
    site_indices = CartesianIndices(Tuple(dims))
    Mf_size = size(cbc.coeff_tensor, N + 1)
    translated_atoms = Vector{Int}(undef, N)
    atoms_sorted_buf = Vector{Int}(undef, N)
    searched_pairs = Set{Vector{Int}}()

    @inbounds for itrans in symmetry.symnum_translation
        atom_site_idx = 0
        for site_idx in 1:N
            ta = symmetry.map_sym[cbc.atoms[site_idx], itrans]
            translated_atoms[site_idx] = ta
            if ta == atom; atom_site_idx = site_idx; end
        end
        copyto!(atoms_sorted_buf, translated_atoms)
        sort!(atoms_sorted_buf)
        atoms_sorted = copy(atoms_sorted_buf)
        atoms_sorted in searched_pairs && continue
        push!(searched_pairs, atoms_sorted)
        atom_site_idx == 0 && continue

        sh_values = Vector{Vector{Float64}}(undef, N)
        l_atom = cbc.ls[atom_site_idx]
        atom_grad_values = Vector{SVector{3,Float64}}(undef, 2*l_atom + 1)

        for (site_idx, translated_atom) in enumerate(translated_atoms)
            l = cbc.ls[site_idx]
            uvec = SVector{3,Float64}(spin_directions[:, translated_atom]...)
            if site_idx == atom_site_idx
                Y, dY = compute_with_gradients(sph, uvec)
                sh_values[site_idx] = [Y[_sph_idx(l, mi)] for mi in 1:(2l+1)]
                for mi in 1:(2l+1)
                    atom_grad_values[mi] = dY[_sph_idx(l, mi)]
                end
            else
                Y = compute(sph, uvec)
                sh_values[site_idx] = [Y[_sph_idx(l, mi)] for mi in 1:(2l+1)]
            end
        end

        grad_result = MVector{3,Float64}(0.0, 0.0, 0.0)
        for mf_idx in 1:Mf_size
            mf_grad = MVector{3,Float64}(0.0, 0.0, 0.0)
            for site_idx_tuple in site_indices
                product = 1.0
                for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
                    site_idx == atom_site_idx && continue
                    product *= sh_values[site_idx][m_idx]
                end
                m_idx_atom = site_idx_tuple.I[atom_site_idx]
                coeff_val = cbc.coeff_tensor[site_idx_tuple.I..., mf_idx]
                mf_grad .+= coeff_val * product .* atom_grad_values[m_idx_atom]
            end
            grad_result .+= cbc.coefficient[mf_idx] .* mf_grad
        end
        result .+= grad_result .* cbc.multiplicity
    end
    return Vector{Float64}(result)
end

# ── sphericart_v2: precompute all harmonics/gradients for all atoms once ───────

function precompute_harmonics_and_grads(
    spin_directions::AbstractMatrix{<:Real},
    sph::SphericalHarmonics,
)
    n_atoms = size(spin_directions, 2)
    n_harm  = (_sph_lmax(sph) + 1)^2
    Y_cache  = Matrix{Float64}(undef, n_atoms, n_harm)
    dY_cache = Matrix{SVector{3,Float64}}(undef, n_atoms, n_harm)
    for atom in 1:n_atoms
        uvec = SVector{3,Float64}(spin_directions[:, atom]...)
        Y, dY = compute_with_gradients(sph, uvec)
        for i in 1:n_harm
            Y_cache[atom, i]  = Y[i]
            dY_cache[atom, i] = dY[i]
        end
    end
    return Y_cache, dY_cache
end

function design_matrix_energy_element_v2(
    cbc,
    symmetry,
    Y_cache::Matrix{Float64},
)::Float64
    result = 0.0
    N = length(cbc.atoms)
    dims = [2*l + 1 for l in cbc.ls]
    site_indices = CartesianIndices(Tuple(dims))
    Mf_size = size(cbc.coeff_tensor, N + 1)
    searched_pairs = Set{Vector{Int}}()

    for itrans in symmetry.symnum_translation
        translated_atoms = [symmetry.map_sym[atom, itrans] for atom in cbc.atoms]
        atoms_sorted = sort(translated_atoms)
        atoms_sorted in searched_pairs && continue
        push!(searched_pairs, atoms_sorted)

        sh_values = [
            [Y_cache[translated_atoms[si], _sph_idx(cbc.ls[si], mi)]
             for mi in 1:(2*cbc.ls[si]+1)]
            for si in 1:N
        ]

        tensor_result = 0.0
        for mf_idx in 1:Mf_size
            mf_contribution = 0.0
            for site_idx_tuple in site_indices
                product = 1.0
                for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
                    product *= sh_values[site_idx][m_idx]
                end
                mf_contribution += cbc.coeff_tensor[site_idx_tuple.I..., mf_idx] * product
            end
            tensor_result += cbc.coefficient[mf_idx] * mf_contribution
        end
        result += tensor_result * cbc.multiplicity
    end
    return result
end

function calc_∇ₑu_v2(
    cbc,
    atom::Integer,
    symmetry,
    Y_cache::Matrix{Float64},
    dY_cache::Matrix{SVector{3,Float64}},
)::Vector{Float64}
    result = MVector{3,Float64}(0.0, 0.0, 0.0)
    N = length(cbc.atoms)
    dims = [2*l + 1 for l in cbc.ls]
    site_indices = CartesianIndices(Tuple(dims))
    Mf_size = size(cbc.coeff_tensor, N + 1)
    translated_atoms = Vector{Int}(undef, N)
    atoms_sorted_buf = Vector{Int}(undef, N)
    searched_pairs = Set{Vector{Int}}()

    @inbounds for itrans in symmetry.symnum_translation
        atom_site_idx = 0
        for site_idx in 1:N
            ta = symmetry.map_sym[cbc.atoms[site_idx], itrans]
            translated_atoms[site_idx] = ta
            if ta == atom; atom_site_idx = site_idx; end
        end
        copyto!(atoms_sorted_buf, translated_atoms)
        sort!(atoms_sorted_buf)
        atoms_sorted = copy(atoms_sorted_buf)
        atoms_sorted in searched_pairs && continue
        push!(searched_pairs, atoms_sorted)
        atom_site_idx == 0 && continue

        sh_values = [
            [Y_cache[translated_atoms[si], _sph_idx(cbc.ls[si], mi)]
             for mi in 1:(2*cbc.ls[si]+1)]
            for si in 1:N
        ]
        l_atom = cbc.ls[atom_site_idx]
        atom_dY = [dY_cache[atom, _sph_idx(l_atom, mi)] for mi in 1:(2*l_atom+1)]

        grad_result = MVector{3,Float64}(0.0, 0.0, 0.0)
        for mf_idx in 1:Mf_size
            mf_grad = MVector{3,Float64}(0.0, 0.0, 0.0)
            for site_idx_tuple in site_indices
                product = 1.0
                for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
                    site_idx == atom_site_idx && continue
                    product *= sh_values[site_idx][m_idx]
                end
                m_idx_atom = site_idx_tuple.I[atom_site_idx]
                coeff_val = cbc.coeff_tensor[site_idx_tuple.I..., mf_idx]
                mf_grad .+= coeff_val * product .* atom_dY[m_idx_atom]
            end
            grad_result .+= cbc.coefficient[mf_idx] .* mf_grad
        end
        result .+= grad_result .* cbc.multiplicity
    end
    return Vector{Float64}(result)
end

# ── full torque matrix with v2 precomputed cache ───────────────────────────────

function build_design_matrix_torque_v2(
    salc_list,
    spinconfig_list,
    num_atoms::Integer,
    symmetry,
    sph::SphericalHarmonics,
)::Matrix{Float64}
    num_salcs       = length(salc_list)
    num_spinconfigs = length(spinconfig_list)
    design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)

    for sc_idx in 1:num_spinconfigs
        spinconfig = spinconfig_list[sc_idx]
        Y_cache, dY_cache = precompute_harmonics_and_grads(spinconfig.spin_directions, sph)
        torque_design_block = zeros(Float64, 3*num_atoms, num_salcs)
        for iatom in 1:num_atoms
            @views dir_iatom = spinconfig.spin_directions[:, iatom]
            for (salc_idx, key_group) in enumerate(salc_list)
                group_grad = MVector{3,Float64}(0.0, 0.0, 0.0)
                for cbc in key_group
                    grad_u = calc_∇ₑu_v2(cbc, iatom, symmetry, Y_cache, dY_cache)
                    group_grad .+= grad_u
                end
                n_C = length(key_group[1].atoms)
                scaling_factor = (4π)^(n_C / 2)
                @views torque_design_block[(3*(iatom-1)+1):(3*iatom), salc_idx] =
                    cross(dir_iatom, Vector{Float64}(group_grad)) * scaling_factor
            end
        end
        design_matrix_list[sc_idx] = torque_design_block
    end
    return vcat(design_matrix_list...)
end

# ── display helper ─────────────────────────────────────────────────────────────

function print_comparison(label::String, b_base, b_sph)
    t_base = median(b_base).time
    t_sph  = median(b_sph).time
    ratio  = t_base / t_sph
    winner = ratio >= 1 ? "SpheriCart          " : "MySphericalHarmonics"
    println("  [$label]")
    @printf "    baseline (MySphericalHarmonics) : %8.2f μs\n" t_base / 1e3
    @printf "    SpheriCart                      : %8.2f μs\n" t_sph  / 1e3
    @printf "    → %s is ~%.2f× faster\n\n" winner max(ratio, 1/ratio)
end

# ── main ───────────────────────────────────────────────────────────────────────

function run_bench()
    ns = BENCH_SAMPLES
    ne = BENCH_EVALS

    example_dir = joinpath(@__DIR__, "examples", BENCH_EXAMPLE)
    input_path  = joinpath(example_dir, "input.toml")
    input = TOML.parsefile(input_path)
    # Resolve datafile relative to example_dir
    datafile_raw = input["regression"]["datafile"]
    if !isabspath(datafile_raw)
        input["regression"]["datafile"] = joinpath(example_dir, datafile_raw)
    end

    print("Loading system ... ")
    system = cd(example_dir) do
        System(input, verbosity = false)
    end
    println("done.")

    spinconfig_list = cd(example_dir) do
        Magesty.SpinConfigs.read_embset(input["regression"]["datafile"])
    end

    salc_list  = system.basisset.salc_list
    symmetry   = system.symmetry
    num_atoms  = system.structure.supercell.num_atoms
    sc1        = spinconfig_list[1]

    # Representative CoupledBasis for microbenchmarks (pick from largest key group)
    cbc_sample = salc_list[end][1]
    iatom_sample = cbc_sample.atoms[1]

    lmax_global = maximum(
        maximum(cbc.ls) for key_group in salc_list for cbc in key_group
    )
    sph = SphericalHarmonics(lmax_global)

    # Precomputed cache for v2 microbenchmarks
    Y_pre, dY_pre = precompute_harmonics_and_grads(sc1.spin_directions, sph)

    println("=" ^ 64)
    println("  Optimize.jl: MySphericalHarmonics vs SpheriCart")
    @printf "  system: %d atoms, lmax = %d, %d SALCs, %d spinconfigs\n" num_atoms lmax_global length(salc_list) length(spinconfig_list)
    @printf "  samples = %d, evals = %d\n" ns ne
    println("=" ^ 64)

    # warm-up
    print("Warming up ... ")
    design_matrix_energy_element(cbc_sample, sc1.spin_directions, symmetry)
    design_matrix_energy_element_v1(cbc_sample, sc1.spin_directions, symmetry, sph)
    design_matrix_energy_element_v2(cbc_sample, symmetry, Y_pre)
    calc_∇ₑu(cbc_sample, iatom_sample, sc1.spin_directions, symmetry)
    calc_∇ₑu_v1(cbc_sample, iatom_sample, sc1.spin_directions, symmetry, sph)
    calc_∇ₑu_v2(cbc_sample, iatom_sample, symmetry, Y_pre, dY_pre)
    build_design_matrix_torque(salc_list, [sc1], num_atoms, symmetry)
    build_design_matrix_torque_v2(salc_list, [sc1], num_atoms, symmetry, sph)
    println("done.\n")

    # ── microbenchmark: design_matrix_energy_element ─────────────────────────
    println("─── design_matrix_energy_element (1 CoupledBasis) ───")
    b_e_base = @benchmark design_matrix_energy_element($cbc_sample, $(sc1.spin_directions), $symmetry) samples=ns evals=ne
    b_e_v1   = @benchmark design_matrix_energy_element_v1($cbc_sample, $(sc1.spin_directions), $symmetry, $sph) samples=ns evals=ne
    b_e_v2   = @benchmark design_matrix_energy_element_v2($cbc_sample, $symmetry, $Y_pre) samples=ns evals=ne

    println("  baseline:"); show(stdout, MIME"text/plain"(), b_e_base); println()
    println("  v1 (compute per atom):"); show(stdout, MIME"text/plain"(), b_e_v1); println()
    println("  v2 (precomputed cache):"); show(stdout, MIME"text/plain"(), b_e_v2); println()
    print_comparison("energy element, baseline vs v1", b_e_base, b_e_v1)
    print_comparison("energy element, baseline vs v2", b_e_base, b_e_v2)

    # ── microbenchmark: calc_∇ₑu ─────────────────────────────────────────────
    println("─── calc_∇ₑu (1 CoupledBasis, 1 atom) ───")
    b_g_base = @benchmark calc_∇ₑu($cbc_sample, $iatom_sample, $(sc1.spin_directions), $symmetry) samples=ns evals=ne
    b_g_v1   = @benchmark calc_∇ₑu_v1($cbc_sample, $iatom_sample, $(sc1.spin_directions), $symmetry, $sph) samples=ns evals=ne
    b_g_v2   = @benchmark calc_∇ₑu_v2($cbc_sample, $iatom_sample, $symmetry, $Y_pre, $dY_pre) samples=ns evals=ne

    println("  baseline:"); show(stdout, MIME"text/plain"(), b_g_base); println()
    println("  v1 (compute_with_gradients per atom):"); show(stdout, MIME"text/plain"(), b_g_v1); println()
    println("  v2 (precomputed cache):"); show(stdout, MIME"text/plain"(), b_g_v2); println()
    print_comparison("calc_∇ₑu, baseline vs v1", b_g_base, b_g_v1)
    print_comparison("calc_∇ₑu, baseline vs v2", b_g_base, b_g_v2)

    # ── macrobenchmark: build_design_matrix_torque (1 spinconfig) ────────────
    println("─── build_design_matrix_torque (1 spinconfig) ───")
    b_t_base = @benchmark build_design_matrix_torque($salc_list, [$(sc1)], $num_atoms, $symmetry) samples=ns evals=ne
    b_t_v2   = @benchmark build_design_matrix_torque_v2($salc_list, [$(sc1)], $num_atoms, $symmetry, $sph) samples=ns evals=ne

    println("  baseline:"); show(stdout, MIME"text/plain"(), b_t_base); println()
    println("  v2 (precomputed cache):"); show(stdout, MIME"text/plain"(), b_t_v2); println()
    print_comparison("torque matrix (1 config), baseline vs v2", b_t_base, b_t_v2)

    println("=" ^ 64)
    println("NOTE: v1 = SpheriCart per atom (minimal code change)")
    println("      v2 = SpheriCart precomputed for all atoms once per spinconfig")
    println("           (requires restructuring build_design_matrix_*)")
    println("=" ^ 64)
end

run_bench()
