#!/usr/bin/env julia

using ArgParse
using DelimitedFiles
using LinearAlgebra
using Interpolations
using Plots

"""
Parse command-line arguments using ArgParse.
"""
function parse_command_line()
    s = ArgParseSettings()
    @add_arg_table s begin
        "file"
            help = "CHGCAR file path"
            arg_type = String
            required = true
        "atom1"
            help = "Index of the first atom"
            arg_type = Int
            required = true
        "atom2"
            help = "Index of the second atom"
            arg_type = Int
            required = true
        "--npoints", "-n"
            help = "Number of sampling points"
            arg_type = Int
            default = 200
        "--ymax"
            help = "Maximum value for y-axis"
            arg_type = Float64
            default = nothing
    end
    return parse_args(s)
end

"""
Read a VASP CHGCAR (or similar) file and extract lattice, atom fractional coordinates,
and the 3D electron density grid.
"""
function read_CHGCAR(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    open(filename) do io
        try
            # Read comment line
            readline(io)
            
            # Read scale factor
            scale_line = strip(readline(io))
            if isempty(scale_line)
                error("Empty scale factor line")
            end
            scale = parse(Float64, scale_line)
            
            # Read lattice vectors
            lat = zeros(3, 3)
            for i in 1:3
                lat_line = strip(readline(io))
                if isempty(lat_line)
                    error("Empty lattice vector line $i")
                end
                lat[i, :] = parse.(Float64, split(lat_line))
            end
            
            # Skip species names
            readline(io)
            
            # Read atom counts
            counts_line = strip(readline(io))
            if isempty(counts_line)
                error("Empty atom counts line")
            end
            counts = parse.(Int, split(counts_line))
            nAtoms = sum(counts)
            if nAtoms <= 0
                error("Invalid number of atoms: $nAtoms")
            end
            
            # Skip coordinate type
            readline(io)
            
            # Read atom fractional coordinates
            atom_frac = zeros(nAtoms, 3)
            for i in 1:nAtoms
                atom_line = strip(readline(io))
                if isempty(atom_line)
                    error("Empty atom coordinate line $i")
                end
                coords = parse.(Float64, split(atom_line))
                if length(coords) < 3
                    error("Insufficient coordinates for atom $i")
                end
                atom_frac[i, :] = coords[1:3]
            end
            
            # Read grid dimensions (skip empty lines)
            dims_line = ""
            while isempty(dims_line) && !eof(io)
                dims_line = strip(readline(io))
            end
            if isempty(dims_line)
                error("Could not find grid dimensions line")
            end
            dims = parse.(Int, split(dims_line))
            if length(dims) != 3
                error("Invalid grid dimensions: expected 3 values, got $(length(dims))")
            end
            nx, ny, nz = dims
            if any(x -> x <= 0, dims)
                error("Invalid grid dimensions: $dims")
            end
            
            # Read density data (skip empty lines)
            ngrid = nx * ny * nz
            data = Float64[]
            while length(data) < ngrid && !eof(io)
                line = strip(readline(io))
                if !isempty(line)
                    try
                        append!(data, parse.(Float64, split(line)))
                    catch e
                        @warn "Skipping invalid line: $line"
                        continue
                    end
                end
            end
            
            if length(data) != ngrid
                error("Grid size mismatch: expected $ngrid, got $(length(data))")
            end
            
            density = reshape(data, nx, ny, nz)
            return scale .* lat, atom_frac, density, (nx, ny, nz)
            
        catch e
            if isa(e, ArgumentError)
                error("Invalid data format in file: $filename - $(e.msg)")
            else
                rethrow(e)
            end
        end
    end
end

"""
Sample density along the line between atoms i and j.
Returns distances (Å) and corresponding densities.
"""
function get_line_density(lat::Array{Float64,2}, atom_frac::Array{Float64,2},
                          density::Array{Float64,3}, dims::NTuple{3,Int},
                          i::Int, j::Int; npoints::Int=200)
    r1 = lat * atom_frac[i, :]
    r2 = lat * atom_frac[j, :]
    ts = range(0.0, 1.0, length=npoints)
    invlat = inv(lat)
    nx, ny, nz = dims
    itp = interpolate(density, BSpline(Linear()), OnGrid())
    gitp = extrapolate(itp, Periodic())
    distances = zeros(npoints)
    dens_vals = zeros(npoints)
    for (k, t) in enumerate(ts)
        pos = (1 - t) * r1 + t * r2
        frac = invlat * pos
        frac .-= floor.(frac)
        x = 1 + frac[1] * (nx - 1)
        y = 1 + frac[2] * (ny - 1)
        z = 1 + frac[3] * (nz - 1)
        dens_vals[k] = gitp(x, y, z)
        distances[k] = norm(pos - r1)
    end
    return distances, dens_vals
end

"""
Main entry point: parse args, compute density profile, and display a plot.
"""
function main()
    args = parse_command_line()
    filename = args["file"]
    i = args["atom1"]
    j = args["atom2"]
    npts = args["npoints"]
    ymax = args["ymax"]

    lattice, atom_frac, density, dims = read_CHGCAR(filename)
    dists, dens = get_line_density(lattice, atom_frac, density, dims, i, j; npoints=npts)

    # Calculate and output total distance
    total_distance = dists[end]
    println("Total distance between Atom $i and Atom $j: $(round(total_distance, digits=5)) Å")

    plt = plot(dists, dens,
               xlabel = "Distance (Å)",
               ylabel = "Electron density",
               title = "Electron density along atom $i → atom $j",
               legend = false,
               xlims = (0.0, total_distance))
    
    # Set y-axis limits
    if ymax !== nothing
        # Set y-axis range: minimum = 0, maximum = specified value
        ylims!(plt, 0.0, ymax)
    else
        # Set y-axis minimum to 0, maximum auto
        current_ylims = ylims(plt)
        ylims!(plt, 0.0, current_ylims[2])
    end
    
    gui()
    readline()
end

# Run script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end