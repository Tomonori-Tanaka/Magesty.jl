using TOML
using LinearAlgebra
using Plots
include("../src/Magesty.jl")


function read_jij(jij_file::AbstractString)
    # Read Jij values from file
    # Format: distance Jij
    # Returns: Vector{Vector{Float64}} with [distance, Jij] pairs
    
    if !isfile(jij_file)
        error("Jij file not found: $jij_file")
    end
    
    data = Vector{Vector{Float64}}()
    
    open(jij_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line) && !startswith(line, "#")
                values = parse.(Float64, split(line))
                if length(values) == 2
                    push!(data, values)
                else
                    error("Invalid line format: $line")
                end
            end
        end
    end
    
    if isempty(data)
        error("No valid data found in file: $jij_file")
    end
    
    println("Jij values:")
    for (i, pair) in enumerate(data)
        println("$(i): distance=$(pair[1]), Jij=$(pair[2])")
    end
    println("--------------------------------")
    
    return data
end

function read_struct(struct_file::AbstractString)
    # Read lattice vectors from file
    # Format: 3x3 matrix of lattice vectors
    # Returns: Matrix{Float64} with 3x3 lattice vectors
    
    if !isfile(struct_file)
        error("Structure file not found: $struct_file")
    end
    
    lattice_vectors = zeros(Float64, 3, 3)
    row_count = 0
    
    open(struct_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line) && !startswith(line, "#")
                values = parse.(Float64, split(line))
                if length(values) == 3
                    row_count += 1
                    if row_count <= 3
                        lattice_vectors[row_count, :] = values
                    else
                        error("Too many rows in lattice vectors file")
                    end
                else
                    error("Invalid line format: $line")
                end
            end
        end
    end
    
    if row_count != 3
        error("Expected 3 rows for lattice vectors, got $row_count")
    end
    
    println("Lattice vectors:")
    println(lattice_vectors)
    println("--------------------------------")
    
    return lattice_vectors
end

function read_kpath(kpath_file::AbstractString)
    # Read k-path from file
    # Format: two k-points (start and end)
    # Returns: Vector of k-points interpolated between start and end
    
    if !isfile(kpath_file)
        error("K-path file not found: $kpath_file")
    end
    
    k_points = Vector{Vector{Float64}}()
    
    open(kpath_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line) && !startswith(line, "#")
                values = parse.(Float64, split(line))
                if length(values) == 3
                    push!(k_points, values)
                else
                    error("Invalid line format: $line")
                end
            end
        end
    end
    
    if length(k_points) != 2
        error("Expected 2 k-points, got $(length(k_points))")
    end
    
    # Interpolate between the two k-points with 50 divisions
    k_start = k_points[1]
    k_end = k_points[2]
    num_divisions = 50
    
    k_path = Vector{Vector{Float64}}()
    for i in 0:num_divisions
        t = i / num_divisions
        k_interpolated = k_start .+ t .* (k_end .- k_start)
        push!(k_path, k_interpolated)
    end
    
    println("K-path:")
    println("Start: $(k_start)")
    println("End: $(k_end)")
    println("Number of points: $(length(k_path))")
    println("--------------------------------")
    
    return k_path
end

function main(input_file::AbstractString, jij_file::AbstractString, struct_file::AbstractString, kpath_file::AbstractString)
    input = TOML.parse(
        open(input_file, "r"),
    )

    config::Magesty.Config4System = Magesty.Config4System(input)
    structure::Magesty.Structure = Magesty.Structure(config, verbosity = true)
    symmetry::Magesty.Symmetry = Magesty.Symmetry(structure, config, verbosity = true)
    cluster::Magesty.Cluster = Magesty.Cluster(structure, symmetry, config, verbosity = true)

    jij::Vector{Vector{Float64}} = read_jij(jij_file)
    lattice_vectors::Matrix{Float64} = read_struct(struct_file)
    reciprocal_lattice_vectors::Matrix{Float64} = 2π * inv(lattice_vectors)
    kpath::Vector{Vector{Float64}} = read_kpath(kpath_file)


    num_atoms::Int = structure.supercell.num_atoms
    center_atom::Int = 1
    j_q_list = zeros(ComplexF64, length(kpath))
    for (k_idx, k) in enumerate(kpath)
        k_cart = reciprocal_lattice_vectors * k
        println("K-point $(k_idx): $(k_cart)")

        j_q = 0.0 + 0.0im
        for i in 1:num_atoms
            if i == center_atom
                continue
            end
            dist_info_list = cluster.min_distance_pairs[center_atom, i]
            for dist_info in dist_info_list
                found = false
                for jij_pair in jij
                    if isapprox(dist_info.distance, jij_pair[1], atol=1e-4)
                        j_q += jij_pair[2] * exp(1.0im * dot(k_cart, dist_info.relative_vector))
                        found = true
                        break
                    end
                end
                if !found
                    error("Jij value not found for distance $(dist_info.distance)")
                end
            end
        end
        j_q_list[k_idx] = j_q
    end
    
    # Check imaginary parts and convert to real if they are sufficiently small
    tolerance::Float64 = 1e-10
    for i in eachindex(j_q_list)
        if abs(imag(j_q_list[i])) < tolerance
            j_q_list[i] = real(j_q_list[i])
        end
    end
    
    println("J_q values (after imaginary part check):")
    for (i, j_q) in enumerate(j_q_list)
        println("$(i): $(j_q)")
    end
    
    # Calculate spin wave dispersion w(k) = 4*(J_0 - J_q)
    J_0 = j_q_list[1]  # J_q at k = 0
    w_k = 2.0 .* (J_0 .- j_q_list)
    
    # Convert w_k to real numbers for plotting
    w_k_real = real.(w_k)
    
    println("Spin wave dispersion w(k):")
    for (i, w) in enumerate(w_k_real)
        println("$(i): w = $(w)")
    end
    
    # Create k-point labels for x-axis
    k_labels = ["Γ" => 1, "X" => length(kpath)]
    
    # Plot spin wave dispersion
    p = plot(1:length(kpath), w_k_real, 
             xlabel="k-point", 
             ylabel="w(k) = 4(J₀ - J_q) [meV]", 
             title="Spin Wave Dispersion",
             marker=:circle,
             linewidth=2,
             legend=false,
             grid=true)
    
    # Add k-point labels
    xticks!(p, [1, length(kpath)], ["Γ", "H"])
    
    # Save plot
    savefig(p, "spinwave_dispersion.png")
    println("Plot saved as spinwave_dispersion.png")
    
    # Display plot
    display(p)
    println("Press Enter to exit the program.")
    readline()
end

input_file = ARGS[1]
jij_file = ARGS[2]
struct_file = ARGS[3]
kpath_file = ARGS[4]

main(input_file, jij_file, struct_file, kpath_file)