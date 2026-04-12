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

function collect_jij(
	input::String,
	reference_atom::Int,
)::Vector{Tuple{Float64, Float64, String}}

	t_setup = @elapsed begin
		structure = Magesty.Structure(input, verbosity = false)
		symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)
		pre = ExchangeTensor.precompute_xml(input)
	end
	println(@sprintf("  [timing] setup (Structure + Symmetry + XML precompute): %.3f s", t_setup))

	num_atoms = structure.supercell.num_atoms

	# Element symbol of the reference atom (i)
	ref_element_index = structure.supercell.kd_int_list[reference_atom]
	ref_element_symbol = structure.kd_name[ref_element_index]

	# Array to store (distance, Jij, element) tuples
	distance_jij_pairs = Vector{Tuple{Float64, Float64, String}}()

	t_loop = @elapsed begin
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
				distance = norm(
					structure.x_image_cart[:, reference_atom, 1] -
					structure.x_image_cart[:, j_atom, cell_index]
				)
				if distance < min_distance
					min_distance = distance
				end
			end

			# Calculate Jij using precomputed data (no XML re-parse)
			jij = ExchangeTensor.convert2tensor_fast(symmetry, reference_atom, j_atom, pre)

			# Get element of the partner atom (j_atom)
			partner_index = structure.supercell.kd_int_list[j_atom]
			partner_symbol = structure.kd_name[partner_index]

			pair_label = string(ref_element_symbol, "–", partner_symbol)
			push!(distance_jij_pairs, (min_distance, jij.isotropic_jij, pair_label))
		end
	end
	println(@sprintf("  [timing] pair loop (%d pairs): %.3f s", num_atoms - 1, t_loop))
	println(@sprintf("  [timing] total collect_jij: %.3f s", t_setup + t_loop))

	# Sort by distance
	sort!(distance_jij_pairs, by = x -> x[1])

	return distance_jij_pairs
end


function plot_jij(
	distance_jij_pairs::Vector{Tuple{Float64, Float64, String}};
	invert_sign::Bool = false,
	half_jij::Bool = false,
	ymin::Real = NaN,
	ymax::Real = NaN,
	markersize::Int = 5,
	no_legend::Bool = false,
)
	if isempty(distance_jij_pairs)
		println("No data to plot.")
		return
	end
	
	# Separate distances and Jij values
	distances = [pair[1] for pair in distance_jij_pairs]
	jij_values = [pair[2] for pair in distance_jij_pairs]
	elements = [pair[3] for pair in distance_jij_pairs]
	
	# Invert sign if requested
	if invert_sign
		jij_values = -jij_values
	end

	# Optionally scale Jij by 1/2
	if half_jij
		jij_values = jij_values ./ 2
	end
	
	# Calculate default y-axis limits rounded to multiples of 10
	jij_min = minimum(jij_values)
	jij_max = maximum(jij_values)

	auto_ymin = floor(jij_min / 10) * 10
	auto_ymax = ceil(jij_max / 10) * 10

	# If user did not specify ymin/ymax, ensure range is at least 20
	if !isfinite(ymin) && !isfinite(ymax) && auto_ymax - auto_ymin < 20
		center = (jij_min + jij_max) / 2
		auto_ymin = floor(center / 10) * 10 - 10
		auto_ymax = ceil(center / 10) * 10 + 10
	end

	# Decide final y-limits (user-specified values override auto limits if finite)
	final_ymin = isfinite(ymin) ? float(ymin) : auto_ymin
	final_ymax = isfinite(ymax) ? float(ymax) : auto_ymax
	
	# Print data to stdout
	println("\nDistance (Å)    Jij (meV)   Pair(i–j)")
	println("------------------------------------------------")
	for (dist, jij, elem) in distance_jij_pairs
		# Apply sign inversion for display if needed
		display_jij = invert_sign ? -jij : jij
		# Apply half scaling if requested
		display_jij = half_jij ? display_jij / 2 : display_jij
		println(@sprintf("%12.6f  %12.6f  %s", dist, display_jij, elem))
	end
	println("------------------------------------------------")
	println(@sprintf("Total pairs: %d", length(distance_jij_pairs)))
	println()
	
	# Create scatter plot (color-coded by element)
	p = scatter(
		distances,
		jij_values,
		group = elements,
		markersize = markersize,
		markeralpha = 0.8,
		title = "Isotropic Jij vs Distance",
		xlabel = "Distance (Å)",
		ylabel = "Jij (meV)",
		xtickfont = font(13),
		ytickfont = font(13),
		legend = no_legend ? false : :topright,
		legendfontsize = 13,
		grid = true,
		size = (800, 600),
		dpi = 300,
		framestyle = :box,
		ylims = (final_ymin, final_ymax)
	)
	
	# Add horizontal line at y=0 (no legend entry)
	hline!(p, [0.0], color = :black, linestyle = :dash, linewidth = 1, alpha = 0.5, label = "")
	
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

        "--half-jij", "-H"
        help = "Plot Jij/2 instead of Jij"
        action = :store_true

        "--ymin"
        help = "Minimum value of y-axis"
        arg_type = Float64
        required = false
        default = NaN

        "--ymax"
        help = "Maximum value of y-axis"
        arg_type = Float64
        required = false
        default = NaN

        "--markersize", "-m"
        help = "Marker size (default: 5)"
        arg_type = Int
        default = 5

        "--no-legend"
        help = "Hide the legend"
        action = :store_true
    end
    args = parse_args(ARGS, s)
    distance_jij_pairs = collect_jij(args["input"], args["reference_atom"])
    plot_jij(
        distance_jij_pairs;
        invert_sign = args["invert-sign"],
        half_jij = args["half-jij"],
        ymin = args["ymin"],
        ymax = args["ymax"],
        markersize = args["markersize"],
        no_legend = args["no-legend"],
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end