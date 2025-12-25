"""
    jij_plot.jl

This script plots isotropic Jij vs distance for a given system.
"""

using LinearAlgebra
using EzXML
using Plots
using ArgParse
using Printf

if !@isdefined(Magesty)
	include("../src/Magesty.jl")
end
using .Magesty

include("./convert2tensor.jl")
using .ExchangeTensor

function collect_jij(input::String, reference_atom::Int)::Vector{Tuple{Float64, Float64}}
	
	structure = Magesty.Structure(input, verbosity = false)
	symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)

	doc = readxml(input)

	num_atoms = structure.supercell.num_atoms
	
	# Array to store distance-Jij pairs
	distance_jij_pairs = Vector{Tuple{Float64, Float64}}()

	for j_atom in 1:num_atoms
		if j_atom == reference_atom
			continue
		end

		# Calculate minimum distance considering periodic boundary conditions
		min_distance = Inf
		for cell_index in 1:27
			if !structure.exist_image[cell_index]
				continue
			end
			# Distance between reference_atom in center cell (cell 1) and j_atom in cell_index
			distance = norm(
				structure.x_image_cart[:, reference_atom, 1] - 
				structure.x_image_cart[:, j_atom, cell_index]
			)
			if distance < min_distance
				min_distance = distance
			end
		end

		# Calculate Jij
		jij = convert2tensor(input, [reference_atom, j_atom])
		
		# Add (distance, Jij) pair
		push!(distance_jij_pairs, (min_distance, jij.isotropic_jij))
	end

	# Sort by distance
	sort!(distance_jij_pairs, by = x -> x[1])
	
	return distance_jij_pairs
end


function plot_jij(distance_jij_pairs::Vector{Tuple{Float64, Float64}}; invert_sign::Bool = false)
	if isempty(distance_jij_pairs)
		println("No data to plot.")
		return
	end
	
	# Separate distances and Jij values
	distances = [pair[1] for pair in distance_jij_pairs]
	jij_values = [pair[2] for pair in distance_jij_pairs]
	
	# Invert sign if requested
	if invert_sign
		jij_values = -jij_values
	end
	
	# Calculate y-axis limits rounded to multiples of 10
	jij_min = minimum(jij_values)
	jij_max = maximum(jij_values)
	
	# Round down to nearest multiple of 10 for minimum
	ymin = floor(jij_min / 10) * 10
	# Round up to nearest multiple of 10 for maximum
	ymax = ceil(jij_max / 10) * 10
	
	# Ensure range is at least 20 if data range is very small
	if ymax - ymin < 20
		center = (jij_min + jij_max) / 2
		ymin = floor(center / 10) * 10 - 10
		ymax = ceil(center / 10) * 10 + 10
	end
	
	# Print data to stdout
	println("\nDistance (Å)    Jij (meV)")
	println("---------------------------")
	for (dist, jij) in distance_jij_pairs
		# Apply sign inversion for display if needed
		display_jij = invert_sign ? -jij : jij
		println(@sprintf("%12.6f  %12.6f", dist, display_jij))
	end
	println("---------------------------")
	println(@sprintf("Total pairs: %d", length(distance_jij_pairs)))
	println()
	
	# Create scatter plot
	p = plot(
		distances,
		jij_values,
		seriestype = :scatter,
		markersize = 5,
		markeralpha = 0.8,
		title = "Isotropic Jij vs Distance",
		xlabel = "Distance (Å)",
		ylabel = "Jij (meV)",
		legend = false,
		grid = true,
		size = (800, 600),
		dpi = 300,
		framestyle = :box,
		ylims = (ymin, ymax)
	)
	
	# Add horizontal line at y=0
	hline!(p, [0.0], color = :black, linestyle = :dash, linewidth = 1, alpha = 0.5)
	
	display(p)
	println("\nPress Enter to exit...")
	readline()
	
	return nothing
end

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "input"
        help = "Input file"
        required = true
        arg_type = String
        
        "reference_atom"
        help = "Reference atom"
        required = true
        arg_type = Int
        
        "--invert-sign", "-i"
        help = "Invert the sign of Jij values"
        action = :store_true
    end
    args = parse_args(ARGS, s)
    distance_jij_pairs = collect_jij(args["input"], args["reference_atom"])
    plot_jij(distance_jij_pairs; invert_sign = args["invert-sign"])
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end