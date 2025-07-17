#!/usr/bin/env julia
module BaderRWIGS

using LinearAlgebra

"""
    ParsedAtIndex

Holds parsed data from a Bader AtIndex.dat file.
"""
struct ParsedAtIndex
    scale::Float64
    lattice::Matrix{Float64}           # 3×3 lattice vectors in Å (transposed for easy conversion)
    atom_positions_frac::Matrix{Float64}  # 3×natoms fractional positions
    AtIndex_grid::Array{Int,3}         # Nx×Ny×Nz grid of atom indices
end

"""
    read_atindex(filename::String) -> ParsedAtIndex

Read a Bader AtIndex.dat file and return parsed data.
"""
function read_atindex(filename::String)::ParsedAtIndex
    open(filename) do io
        readline(io)  # title
        scale = parse(Float64, strip(readline(io)))
        # lattice (read as 3×3, will be transposed for easy conversion)
        lattice = zeros(3,3)
        for i in 1:3
            lattice[i, :] = parse.(Float64, split(strip(readline(io))))
        end
        # Transpose lattice for easy conversion: lattice * atom_positions_frac
        lattice = lattice'
        # skip species line
        split(strip(readline(io)))
        counts = parse.(Int, split(strip(readline(io))))
        natoms = sum(counts)
        _ = strip(readline(io))  # "Direct" or "Cartesian"
        # atom fractional positions (3×natoms matrix)
        atom_positions_frac = zeros(3, natoms)
        for i in 1:natoms
            atom_positions_frac[:, i] = parse.(Float64, split(strip(readline(io))))
        end
        # read grid dims
        line = strip(readline(io))
        while isempty(line)
            line = strip(readline(io))
        end
        dims = parse.(Int, split(line))
        Nx, Ny, Nz = dims
        # read grid values
        flat = Float64[]
        while !eof(io)
            append!(flat, parse.(Float64, split(strip(readline(io)))))
        end
        idx_vals = Int.(round.(flat))
        AtIndex_grid = reshape(idx_vals, Nx, Ny, Nz)
        # convert lattice to Cartesian (already transposed)
        lattice_cart = lattice * scale
        return ParsedAtIndex(scale, lattice_cart, atom_positions_frac, AtIndex_grid)
    end
end

"""
    compute_bader_radii(data::ParsedAtIndex; nsteps::Int=200)

Compute spherical Bader radii under periodic boundary conditions.
Returns a vector of radii (Å) for each atom.
"""
function compute_bader_radii(data::ParsedAtIndex; nsteps::Int=200)
    lat = data.lattice
    pos_frac = data.atom_positions_frac
    grid = data.AtIndex_grid
    Nx, Ny, Nz = size(grid)
    natoms = size(pos_frac, 2)  # 3×natoms matrix
    radii = fill(Inf, natoms)
    
    # Convert all atom positions to Cartesian coordinates at once
    pos_cart = lat * pos_frac  # 3×natoms matrix
    
    # loop atoms
    for iatom in 1:natoms
        p_i = pos_frac[:, iatom]  # 3×1 vector (fractional)
        p_cart_i = pos_cart[:, iatom]  # 3×1 vector (Cartesian)
        for jatom in 1:natoms
            if iatom == jatom
                continue
            end
            p_j = pos_frac[:, jatom]  # 3×1 vector (fractional)
            # minimum-image displacement in fractional coordinates
            δ = p_j .- p_i
            δ .-= round.(δ)  # map into [-0.5,0.5)
            # sample along line
            for step in 1:nsteps
                t = step / nsteps
                p_frac_t = p_i .+ δ .* t
                # wrap into [0,1)
                p_frac_t .-= floor.(p_frac_t)
                # grid indices
                ix = clamp(Int(fld(p_frac_t[1]*Nx,1)) + 1, 1, Nx)
                iy = clamp(Int(fld(p_frac_t[2]*Ny,1)) + 1, 1, Ny)
                iz = clamp(Int(fld(p_frac_t[3]*Nz,1)) + 1, 1, Nz)
                if grid[ix, iy, iz] != iatom
                    # crossing point found
                    p_cart_t = lat * p_frac_t
                    r = norm(p_cart_t - p_cart_i)
                    radii[iatom] = min(radii[iatom], r)
                    break
                end
            end
        end
    end
    return radii
end

end  # module BaderRWIGS

# script entry
if abspath(PROGRAM_FILE) == @__FILE__
    using .BaderRWIGS
    if length(ARGS) < 1
        println("Usage: julia bader_rwigscode.jl <AtIndex.dat> [nsteps]")
        exit(1)
    end
    file = ARGS[1]
    nsteps = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    data = BaderRWIGS.read_atindex(file)
    radii = BaderRWIGS.compute_bader_radii(data; nsteps=nsteps)
    println("Atom radii (Å):")
    for (i, r) in enumerate(radii)
        println("  Atom $(i): $(round(r,digits=4))")
    end
    println("RWIGS = ", join(round.(radii;digits=4)," "))
end