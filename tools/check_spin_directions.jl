#!/usr/bin/env julia

using LinearAlgebra
using Plots
gr()  # Set GR as the backend
using ArgParse

include("vasptools.jl")
using .VaspTools

"""
    plot_spin_directions(spin_matrix::AbstractMatrix{<:Real})

Plot normalized spin vectors in 3D space.
All vectors are normalized and plotted from the origin.

# Arguments
- `spin_matrix`: 3×N matrix containing N spin vectors
"""
function plot_spin_directions(spin_matrix::AbstractMatrix{<:Real})
    # Normalize all vectors
    normalized_vectors = zeros(size(spin_matrix))
    for (i, col) in enumerate(eachcol(spin_matrix))
        if !isapprox(norm(col), 0.0, atol=1e-10)
            normalized_vectors[:, i] = normalize(col)
        end
    end

    # Create 3D plot
    p = plot(
        title="Normalized Spin Directions",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        aspect_ratio=:equal,
        camera=(30, 30),  # Initial camera angle
        size=(800, 800),  # Larger window size
    )

    # Plot unit sphere
    θ = range(0, 2π, length=100)
    φ = range(0, π, length=50)
    x = [sin(φ) * cos(θ) for φ in φ, θ in θ]
    y = [sin(φ) * sin(θ) for φ in φ, θ in θ]
    z = [cos(φ) for φ in φ, θ in θ]
    surface!(p, x, y, z, alpha=0.1, color=:lightblue)

    # Plot vectors
    for (i, col) in enumerate(eachcol(normalized_vectors))
        if !isapprox(norm(col), 0.0, atol=1e-10)
            quiver!(p, 
                [0], [0], [0],  # Start points
                quiver=([col[1]], [col[2]], [col[3]]),  # Vector components
                color=:red,
                label=i == 1 ? "Spin vectors" : ""
            )
        end
    end

    display(p)  # Explicitly display the plot
    return p
end

function main()
    s = ArgParseSettings(description = """
    Visualize normalized spin vectors in 3D space.
    All vectors are normalized and plotted from the origin.
    """)

    @add_arg_table s begin
        "input"
        help = "Input file containing spin vectors (INCAR file)"
        required = true
    end

    args = parse_args(s)
    
    # Read the input file
    incar_dict = parse_incar(args["input"])
    spin_list = incar_dict[:MAGMOM]
    spin_matrix = reshape(spin_list, 3, length(spin_list) ÷ 3)
    
    # Check if the matrix has the correct dimensions
    if size(spin_matrix, 1) != 3
        error("Input matrix must have 3 rows")
    end
    
    # Create and display the plot
    p = plot_spin_directions(spin_matrix)
    
    # Wait for user input before closing
    println("Press Enter to exit...")
    readline()
end

# Execute the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
