using Test
using TOML
using StaticArrays
using Magesty

function reference_calc_grad(
    cbc::Magesty.Basis.CoupledBasis_with_coefficient,
    atom::Integer,
    spin_directions::AbstractMatrix{<:Real},
    symmetry::Magesty.Symmetry,
)::Vector{Float64}
    result = MVector{3, Float64}(0.0, 0.0, 0.0)
    nsites = length(cbc.atoms)
    dims = [2 * l + 1 for l in cbc.ls]
    site_indices = CartesianIndices(Tuple(dims))
    mf_size = size(cbc.coeff_tensor, nsites + 1)

    searched_pairs = Set{Vector{Int}}()
    for itrans in symmetry.symnum_translation
        translated_atoms = [symmetry.map_sym[a, itrans] for a in cbc.atoms]
        atoms_sorted = sort(translated_atoms)
        if atoms_sorted in searched_pairs
            continue
        end
        push!(searched_pairs, atoms_sorted)

        atom_site_idx = findfirst(==(atom), translated_atoms)
        isnothing(atom_site_idx) && continue

        sh_values = Vector{Vector{Float64}}(undef, nsites)
        l_atom = cbc.ls[atom_site_idx]
        atom_grad_values = Vector{Vector{Float64}}(undef, 2 * l_atom + 1)
        for (site_idx, ta) in enumerate(translated_atoms)
            l = cbc.ls[site_idx]
            sh_values[site_idx] = Vector{Float64}(undef, 2 * l + 1)
            for m_idx in 1:(2 * l + 1)
                m = m_idx - l - 1
                sh_values[site_idx][m_idx] = @views Magesty.MySphericalHarmonics.Zₗₘ(
                    l,
                    m,
                    spin_directions[:, ta],
                )
                if site_idx == atom_site_idx
                    atom_grad_values[m_idx] = @views Magesty.MySphericalHarmonics.∂ᵢZlm(
                        l,
                        m,
                        spin_directions[:, ta],
                    )
                end
            end
        end

        grad_result = MVector{3, Float64}(0.0, 0.0, 0.0)
        for mf_idx in 1:mf_size
            mf_grad_contrib = MVector{3, Float64}(0.0, 0.0, 0.0)
            for site_idx_tuple in site_indices
                product = 1.0
                for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
                    if site_idx == atom_site_idx
                        continue
                    end
                    product *= sh_values[site_idx][m_idx]
                end
                m_idx_atom = site_idx_tuple.I[atom_site_idx]
                grad_atom = atom_grad_values[m_idx_atom]
                coeff_val = cbc.coeff_tensor[site_idx_tuple.I..., mf_idx]
                mf_grad_contrib .+= coeff_val * product .* grad_atom
            end
            grad_result .+= cbc.coefficient[mf_idx] .* mf_grad_contrib
        end

        result .+= grad_result .* cbc.multiplicity
    end
    return Vector{Float64}(result)
end

function reference_torque_block(
    salc_list::AbstractVector{Vector{Magesty.Basis.CoupledBasis_with_coefficient}},
    spinconfig::Magesty.SpinConfigs.SpinConfig,
    num_atoms::Integer,
    symmetry::Magesty.Symmetry,
)::Matrix{Float64}
    num_salcs = length(salc_list)
    block = zeros(Float64, 3 * num_atoms, num_salcs)
    for iatom in 1:num_atoms
        dir_iatom = @view spinconfig.spin_directions[:, iatom]
        for (salc_idx, key_group) in enumerate(salc_list)
            group_grad = zeros(3)
            for cbc in key_group
                group_grad .+= reference_calc_grad(cbc, iatom, spinconfig.spin_directions, symmetry)
            end
            n_C = length(key_group[1].atoms)
            scaling = (4 * pi)^(n_C / 2)
            block[(3 * (iatom - 1) + 1):(3 * iatom), salc_idx] .= cross(dir_iatom, group_grad) .* scaling
        end
    end
    return block
end

@testset "Optimize numerical consistency" begin
    input_path = joinpath(@__DIR__, "..", "examples", "fege_2x2x2", "input.toml")
    input = TOML.parse(open(input_path, "r"))
    input["regression"]["datafile"] = joinpath(@__DIR__, "..", "examples", "fege_2x2x2", "EMBSET")

    system = Magesty.System(input, verbosity = false)
    spincluster = Magesty.SpinCluster(system, input, verbosity = false)

    sc = spincluster.optimize.spinconfig_list[1]
    salc_list = spincluster.basisset.salc_list
    symmetry = spincluster.symmetry
    num_atoms = spincluster.structure.supercell.num_atoms

    max_groups = min(length(salc_list), 5)
    for salc_idx in 1:max_groups
        key_group = salc_list[salc_idx]
        max_cbc = min(length(key_group), 2)
        for cbc in key_group[1:max_cbc]
            unique_atoms = unique(cbc.atoms)
            for atom in unique_atoms
                ref = reference_calc_grad(cbc, atom, sc.spin_directions, symmetry)
                now = Magesty.Optimize.calc_∇ₑu(cbc, atom, sc.spin_directions, symmetry)
                @test isapprox(now, ref; atol = 1e-12, rtol = 1e-10)
            end
        end
    end

    ref_torque = reference_torque_block(salc_list, sc, num_atoms, symmetry)
    now_torque = Magesty.Optimize.build_design_matrix_torque(salc_list, [sc], num_atoms, symmetry)
    @test isapprox(now_torque, ref_torque; atol = 1e-12, rtol = 1e-10)
end
