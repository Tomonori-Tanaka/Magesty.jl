#!/usr/bin/env julia
# bader_line_plot.jl
#
# Reads a Bader AtIndex.dat file and plots the atomic index along
# the line connecting two specified atoms under PBC, using ArgParse.

using LinearAlgebra
using DelimitedFiles
using Plots
using ArgParse
using Printf

# --- Parse AtIndex.dat ---
struct ParsedAtIndex
    lattice::Matrix{Float64}    # 3×3 Cartesian lattice (Å)
    frac_pos::Matrix{Float64}   # natoms×3 fractional coordinates
    grid::Array{Int,3}          # AtIndex grid NX×NY×NZ
end

function read_atindex(filename::String)::ParsedAtIndex
    open(filename) do io
        readline(io)  # title
        scale = parse(Float64, strip(readline(io)))
        # read lattice vectors
        lat_frac = [parse.(Float64, split(strip(readline(io)))) for _ in 1:3]
        lattice = hcat(lat_frac...) * scale
        # skip species and counts
        readline(io)
        counts = parse.(Int, split(strip(readline(io))))
        natoms = sum(counts)
        # skip coord type
        readline(io)
        # read fractional positions
        frac_pos = Array{Float64,2}(undef, natoms, 3)
        for i in 1:natoms
            frac_pos[i, :] = parse.(Float64, split(strip(readline(io))))
        end
        # read grid dimensions
        line = strip(readline(io))
        while isempty(line)
            line = strip(readline(io))
        end
        NX, NY, NZ = parse.(Int, split(line))
        # read flat grid values
        vals = Float64[]
        while !eof(io)
            append!(vals, parse.(Float64, split(strip(readline(io)))))
        end
        grid = reshape(Int.(round.(vals)), NX, NY, NZ)
        return ParsedAtIndex(lattice, frac_pos, grid)
    end
end

# Get AtIndex at fractional coordinate (no PBC) - nearest neighbor
grid_index(p_frac, grid) = begin
    NX, NY, NZ = size(grid)
    # No PBC wrapping - use direct fractional coordinates
    ix = clamp(Int(round(p_frac[1]*NX)), 1, NX)
    iy = clamp(Int(round(p_frac[2]*NY)), 1, NY)
    iz = clamp(Int(round(p_frac[3]*NZ)), 1, NZ)
    return grid[ix, iy, iz]
end

# Get AtIndex at fractional coordinate using linear interpolation (no PBC)
function grid_index_interpolated(p_frac::Vector{Float64}, grid::Array{Int,3})::Float64
    NX, NY, NZ = size(grid)
    
    # Convert fractional coordinates to grid coordinates
    x_grid = p_frac[1] * NX
    y_grid = p_frac[2] * NY
    z_grid = p_frac[3] * NZ
    
    # Clamp to grid boundaries
    x_grid = clamp(x_grid, 1.0, Float64(NX))
    y_grid = clamp(y_grid, 1.0, Float64(NY))
    z_grid = clamp(z_grid, 1.0, Float64(NZ))
    
    # Get integer indices for interpolation
    ix_low = Int(floor(x_grid))
    iy_low = Int(floor(y_grid))
    iz_low = Int(floor(z_grid))
    
    ix_high = min(ix_low + 1, NX)
    iy_high = min(iy_low + 1, NY)
    iz_high = min(iz_low + 1, NZ)
    
    # Calculate interpolation weights
    wx = x_grid - ix_low
    wy = y_grid - iy_low
    wz = z_grid - iz_low
    
    # Perform trilinear interpolation
    # Corner values
    v000 = Float64(grid[ix_low, iy_low, iz_low])
    v001 = Float64(grid[ix_low, iy_low, iz_high])
    v010 = Float64(grid[ix_low, iy_high, iz_low])
    v011 = Float64(grid[ix_low, iy_high, iz_high])
    v100 = Float64(grid[ix_high, iy_low, iz_low])
    v101 = Float64(grid[ix_high, iy_low, iz_high])
    v110 = Float64(grid[ix_high, iy_high, iz_low])
    v111 = Float64(grid[ix_high, iy_high, iz_high])
    
    # Interpolate along x-axis
    v00 = v000 * (1 - wx) + v100 * wx
    v01 = v001 * (1 - wx) + v101 * wx
    v10 = v010 * (1 - wx) + v110 * wx
    v11 = v011 * (1 - wx) + v111 * wx
    
    # Interpolate along y-axis
    v0 = v00 * (1 - wy) + v10 * wy
    v1 = v01 * (1 - wy) + v11 * wy
    
    # Interpolate along z-axis
    result = v0 * (1 - wz) + v1 * wz
    
    return result
end

"""
    min_image_disp(p, q)

Calculate the minimum-image displacement between two fractional coordinates under PBC.

# Arguments
- `p::Vector{Float64}`: First fractional coordinate
- `q::Vector{Float64}`: Second fractional coordinate

# Returns
- `Vector{Float64}`: Minimum-image displacement vector

# Notes
- Applies periodic boundary conditions to find the shortest path between points
- Returns displacement vector that minimizes the distance
"""
function min_image_disp(p::Vector{Float64}, q::Vector{Float64})::Vector{Float64}
    δ = q .- p
    # Apply periodic boundary conditions
    δ .-= round.(δ)
    return δ
end

# --- Argument parsing ---
function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "filename"
            help = "Path to the AtIndex.dat file"
        "atom1"
            help = "Index of the first atom (1-based)"
            arg_type = Int
        "atom2"
            help = "Index of the second atom (1-based)"
            arg_type = Int
        "--npoints"
            help = "Number of sampling points along the line"
            arg_type = Int
            default = 200
        "--interpolate"
            help = "Use linear interpolation for AtIndex values"
            action = :store_true
    end
    return parse_args(s)
end

# --- Main ---
function main()
    args = parse_command_line()
    filename = args["filename"]
    atom1 = args["atom1"]
    atom2 = args["atom2"]
    npoints = args["npoints"]
    use_interpolation = args["interpolate"]

    data = read_atindex(filename)
    lattice = data.lattice
    frac = data.frac_pos
    grid = data.grid
    natoms = size(frac, 1)

    # validate indices
    for a in (atom1, atom2)
        if a < 1 || a > natoms
            error("Atom index $a out of range [1,$natoms]")
        end
    end

    # compute direct displacement and distance (no PBC)
    p1_frac = frac[atom1, :]
    p2_frac = frac[atom2, :]
    δ = p2_frac .- p1_frac  # Direct displacement without PBC
    
    # Calculate direct distance
    p1_cart = lattice * p1_frac
    p2_cart = lattice * p2_frac
    total_distance = norm(p2_cart - p1_cart)
    
    println("Total distance between Atom $atom1 and Atom $atom2 (direct): $(round(total_distance, digits=5)) Å")

    # sample along the line
    ts = range(0.0, 1.0, length=npoints)
    dists = zeros(Float64, npoints)
    indices = use_interpolation ? zeros(Float64, npoints) : zeros(Int, npoints)
    for (k, t) in enumerate(ts)
        pf = p1_frac .+ δ .* t
        if use_interpolation
            indices[k] = grid_index_interpolated(pf, grid)
        else
            indices[k] = grid_index(pf, grid)
        end
        # No PBC wrapping - use direct coordinates
        p_cart = lattice * pf
        dists[k] = norm(p_cart - p1_cart)
    end

    # Find points where the index changes
    change_points = Vector{Tuple{Float64, Any, Any}}()  # (distance, old_index, new_index)
    for i in eachindex(indices)[2:end]
        if indices[i] != indices[i-1]
            # Interpolate the exact distance where the change occurs
            t_change = ts[i-1] + (ts[i] - ts[i-1]) * 0.5  # Midpoint approximation
            pf_change = p1_frac .+ δ .* t_change
            # No PBC wrapping - use direct coordinates
            p_change_cart = lattice * pf_change
            dist_change = norm(p_change_cart - p1_cart)
            
            push!(change_points, (dist_change, indices[i-1], indices[i]))
        end
    end

    # Output change points
    if isempty(change_points)
        println("No index changes found along the line.")
    else
        println("\nIndex change points:")
        println("Distance (Å) | Old Index | New Index")
        println("-" ^ 40)
        for (dist, old_idx, new_idx) in change_points
            if use_interpolation
                println(@sprintf("%10.5f | %9.3f | %9.3f", dist, old_idx, new_idx))
            else
                println(@sprintf("%10.5f | %9d | %9d", dist, old_idx, new_idx))
            end
        end
    end

    # plot
    p = plot(dists, indices;
             xlabel = "Distance from Atom $atom1 (Å)",
             ylabel = "AtIndex along line to Atom $atom2",
             title  = "AtIndex profile: Atom $atom1 → Atom $atom2",
             legend = false,
             marker = :circle,
             line = :path)
    
    # Add vertical lines at change points
    for (dist, old_idx, new_idx) in change_points
        vline!([dist], color=:red, style=:dash, alpha=0.7, label="")
    end
    
    gui()
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end