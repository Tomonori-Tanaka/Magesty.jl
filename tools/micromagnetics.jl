"""
Caculation tool for the derivation of the micromagnetics model parameters
"""

using ArgParse
using LinearAlgebra
using EzXML
using JLD2
using Base.Threads

include("../src/Magesty.jl")
using .Magesty
include("./convert2tensor.jl")


function calc_dm_vector(j_ex::Array{Float64, 2})::Vector{Float64}
    x = 1/2*(j_ex[2, 3] - j_ex[3, 2])
    y = 1/2*(j_ex[3, 1] - j_ex[1, 3])
    z = 1/2*(j_ex[1, 2] - j_ex[2, 1])
    return [x, y, z]
end

function calc_micromagnetics(input_xml::String, system::System, cutoff::Union{Float64, Nothing}=nothing)::Tuple{Matrix{Float64}, Matrix{Float64}}
    num_atoms = system.structure.supercell.num_atoms
    atoms_in_prim = system.symmetry.atoms_in_prim   # atom indices in the primitive cell
    min_distance_pairs = system.cluster.min_distance_pairs

    nthreads = Threads.nthreads()
    stiff_loc = [zeros(Float64, 3, 3) for _ in 1:nthreads]
    spiral_loc = [zeros(Float64, 3, 3) for _ in 1:nthreads]

    stiffness_matrix = zeros(Float64, 3, 3)
    spiralization_matrix = zeros(Float64, 3, 3)
    @threads for idx in 1:length(atoms_in_prim)
        i_atom = atoms_in_prim[idx]
        t_id = Threads.threadid()
        stiff  = stiff_loc[t_id]
        spiral = spiral_loc[t_id]

        for i_pair in 1:num_atoms
            if i_pair == i_atom; continue; end

            dist_info_list = min_distance_pairs[i_atom, i_pair]
            for dist_info in dist_info_list
                # Check cutoff distance if specified
                if cutoff !== nothing
                    distance = norm(dist_info.relative_vector)
                    if distance > cutoff
                        continue
                    end
                end

                cell_index = dist_info.cell_index
                # relative vector from i_atom to i_pair in Cartesian coordinates
                relvec::Vector{Float64} = dist_info.relative_vector

                # calculate the exchange interaction tensor
                exchange_tensor = convert2tensor(input_xml, [i_atom, i_pair], cell_index)

                jij = 1/3*tr(exchange_tensor)
                # Avoid creating temporary arrays by using direct element access
                for i in 1:3, j in 1:3
                    stiff[i, j] += 0.5 * jij * relvec[i] * relvec[j]
                end
                
                dm_vector = calc_dm_vector(exchange_tensor)
                # Avoid creating temporary arrays by using direct element access
                for i in 1:3, j in 1:3
                    spiral[i, j] += dm_vector[i] * relvec[j]
                end
            end
        end
        stiff_loc[t_id] = stiff
        spiral_loc[t_id] = spiral
    end
    for t_id in 1:nthreads
        stiffness_matrix .+= stiff_loc[t_id]
        spiralization_matrix .+= spiral_loc[t_id]
    end

    return stiffness_matrix, spiralization_matrix
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_xml", "-x"
        help = "Input xml file"
        required = true

        "--input_jld2", "-j"
        help = "Input jld2 file for System struct"
        required = true

        "--cutoff", "-c"
        help = "Cutoff distance for atom pairs (in Angstrom). If not specified, all atom pairs are considered"
        default = nothing
    end

    args = parse_args(s)
    @load args["input_jld2"] system

    cutoff = args["cutoff"] === nothing ? nothing : parse(Float64, args["cutoff"])
    stiffness_matrix, spiralization_matrix = calc_micromagnetics(args["input_xml"], system, cutoff)
    println("Stiffness matrix:")
    display(stiffness_matrix)
    println("")
    println("Spiralization matrix:")
    display(spiralization_matrix)
end

main()