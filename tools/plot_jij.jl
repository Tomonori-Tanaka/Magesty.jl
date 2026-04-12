"""
    plot_jij.jl

Plot isotropic Jij vs distance for ALL pairs in the system.
Optionally filter by element pair (--element1 and --element2).
Accepts multiple XML files; each file is plotted in a distinct color.
"""

include("./plot_jij_atom.jl")

function collect_jij_all(
	input::String;
	element1::Union{String, Nothing} = nothing,
	element2::Union{String, Nothing} = nothing,
)::Vector{Tuple{Float64, Float64, String}}
	t_setup = @elapsed begin
		structure = Magesty.Structure(input, verbosity = false)
		symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)
		pre = ExchangeTensor.precompute_xml(input)
	end
	println(@sprintf("  [timing] setup (Structure + Symmetry + XML precompute): %.3f s", t_setup))

	num_atoms = structure.supercell.num_atoms

	# Normalize element filter: sorted pair for matching (e.g. "Fe–O" for both Fe,O and O,Fe)
	filter_label = nothing
	if element1 !== nothing && element2 !== nothing
		e1, e2 = strip(element1), strip(element2)
		if isempty(e1) || isempty(e2)
			error("--element1 and --element2 must be non-empty when both specified")
		end
		filter_label = e1 < e2 ? string(e1, "–", e2) : string(e2, "–", e1)
	end

	distance_jij_pairs = Vector{Tuple{Float64, Float64, String}}()

	num_pairs = num_atoms * (num_atoms - 1) ÷ 2
	t_loop = @elapsed begin
	for i in 1:num_atoms
		for j in (i + 1):num_atoms
			# Minimum distance considering periodic boundary conditions
			min_distance = Inf
			for cell_index in 1:27
				if !structure.exist_image[cell_index]
					continue
				end
				distance = norm(
					structure.x_image_cart[:, i, 1] -
					structure.x_image_cart[:, j, cell_index]
				)
				if distance < min_distance
					min_distance = distance
				end
			end

			jij = ExchangeTensor.convert2tensor_fast(symmetry, i, j, pre)

			elem_i = structure.kd_name[structure.supercell.kd_int_list[i]]
			elem_j = structure.kd_name[structure.supercell.kd_int_list[j]]
			# Canonical pair label (alphabetically sorted for consistency)
			pair_label = elem_i < elem_j ? string(elem_i, "–", elem_j) : string(elem_j, "–", elem_i)

			if filter_label !== nothing && pair_label != filter_label
				continue
			end

			push!(distance_jij_pairs, (min_distance, jij.isotropic_jij, pair_label))
		end
	end
	end # t_loop
	println(@sprintf("  [timing] pair loop (%d pairs): %.3f s", num_pairs, t_loop))
	println(@sprintf("  [timing] total collect_jij_all: %.3f s", t_setup + t_loop))

	sort!(distance_jij_pairs, by = x -> x[1])
	return distance_jij_pairs
end

# Palette used when multiple datasets are plotted together.
# Each dataset gets one color; all its data points share that color.
const MULTI_COLORS = [:royalblue, :crimson, :forestgreen, :darkorange, :purple, :saddlebrown, :teal, :deeppink]

"""
Plot multiple datasets on a single figure.  Each dataset is identified by a
label string and rendered in a distinct color.  The y-axis limits are computed
from all data combined.

`datasets` is a vector of `(label, pairs)` where `pairs` is the output of
`collect_jij_all`.
"""
function plot_jij_multi(
	datasets::Vector{Tuple{String, Vector{Tuple{Float64, Float64, String}}}};
	invert_sign::Bool = false,
	half_jij::Bool = false,
	ymin::Real = NaN,
	ymax::Real = NaN,
	markersize::Int = 5,
	no_legend::Bool = false,
)
	all_jij_raw = [pair[2] for (_, pairs) in datasets for pair in pairs]
	if isempty(all_jij_raw)
		println("No data to plot.")
		return
	end

	sign_factor = invert_sign ? -1.0 : 1.0
	half_factor = half_jij ? 0.5 : 1.0
	scale = sign_factor * half_factor

	all_jij = all_jij_raw .* scale
	jij_min = minimum(all_jij)
	jij_max = maximum(all_jij)

	auto_ymin = floor(jij_min / 10) * 10
	auto_ymax = ceil(jij_max / 10) * 10

	if !isfinite(ymin) && !isfinite(ymax) && auto_ymax - auto_ymin < 20
		center = (jij_min + jij_max) / 2
		auto_ymin = floor(center / 10) * 10 - 10
		auto_ymax = ceil(center / 10) * 10 + 10
	end

	final_ymin = isfinite(ymin) ? float(ymin) : auto_ymin
	final_ymax = isfinite(ymax) ? float(ymax) : auto_ymax

	# Print data table
	for (label, pairs) in datasets
		println("\n=== $label ===")
		println("Distance (Å)    Jij (meV)   Pair(i–j)")
		println("------------------------------------------------")
		for (dist, jij, elem) in pairs
			println(@sprintf("%12.6f  %12.6f  %s", dist, jij * scale, elem))
		end
		println("------------------------------------------------")
		println(@sprintf("Total pairs: %d", length(pairs)))
	end
	println()

	p = plot(
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
		ylims = (final_ymin, final_ymax),
	)
	hline!(p, [0.0], color = :black, linestyle = :dash, linewidth = 1, alpha = 0.5, label = "")

	for (idx, (label, pairs)) in enumerate(datasets)
		isempty(pairs) && continue
		distances = [pair[1] for pair in pairs]
		jij_values = [pair[2] for pair in pairs] .* scale
		color = MULTI_COLORS[mod1(idx, length(MULTI_COLORS))]
		scatter!(p, distances, jij_values;
			label = label,
			color = color,
			markersize = markersize,
			markeralpha = 0.8,
		)
	end

	display(p)
	println("\nPress Enter to exit...")
	readline()
	return nothing
end

function main()
	s = ArgParseSettings(
		description = "Plot isotropic Jij vs distance for all pairs. Optionally filter by element pair.",
	)
	@add_arg_table s begin
		"input"
		help = "Input XML file(s)"
		required = true
		arg_type = String
		nargs = '+'

		"--label", "-l"
		help = "Label(s) for each input file (comma-separated, e.g. 'A,B,C'). Defaults to file basename."
		arg_type = String
		default = nothing

		"--element1", "-e"
		help = "First element for filtering (use with --element2)"
		arg_type = String
		default = nothing

		"--element2", "-E"
		help = "Second element for filtering (use with --element1)"
		arg_type = String
		default = nothing

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

	# Validate: if one element is specified, both must be
	e1, e2 = args["element1"], args["element2"]
	if (e1 !== nothing && e2 === nothing) || (e1 === nothing && e2 !== nothing)
		error("--element1 and --element2 must be specified together")
	end

	input_files = args["input"]

	# Build per-file labels
	if args["label"] !== nothing
		labels = strip.(split(args["label"], ","))
		if length(labels) != length(input_files)
			error("Number of --label entries ($(length(labels))) must match number of input files ($(length(input_files)))")
		end
	else
		# Default: use the basename without extension
		labels = [splitext(basename(f))[1] for f in input_files]
	end

	# Collect data per file
	datasets = Vector{Tuple{String, Vector{Tuple{Float64, Float64, String}}}}()
	for (file, label) in zip(input_files, labels)
		println("\n=== Processing: $file (label: $label) ===")
		pairs = collect_jij_all(file; element1 = e1, element2 = e2)
		push!(datasets, (label, pairs))
	end

	if length(datasets) == 1
		# Single file: use the existing plot_jij (color-coded by pair type)
		plot_jij(
			datasets[1][2];
			invert_sign = args["invert-sign"],
			half_jij = args["half-jij"],
			ymin = args["ymin"],
			ymax = args["ymax"],
			markersize = args["markersize"],
			no_legend = args["no-legend"],
		)
	else
		# Multiple files: color-coded by file
		plot_jij_multi(
			datasets;
			invert_sign = args["invert-sign"],
			half_jij = args["half-jij"],
			ymin = args["ymin"],
			ymax = args["ymax"],
			markersize = args["markersize"],
			no_legend = args["no-legend"],
		)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
