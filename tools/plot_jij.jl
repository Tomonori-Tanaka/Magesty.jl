"""
    plot_jij.jl

Plot isotropic Jij vs distance for ALL pairs in the system.
Optionally filter by element pair (--element1 and --element2).
"""

include("./plot_jij_atom.jl")

function collect_jij_all(
	input::String;
	element1::Union{String, Nothing} = nothing,
	element2::Union{String, Nothing} = nothing,
)::Vector{Tuple{Float64, Float64, String}}
	structure = Magesty.Structure(input, verbosity = false)
	symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)

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

			jij = convert2tensor(input, [i, j])

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

	sort!(distance_jij_pairs, by = x -> x[1])
	return distance_jij_pairs
end

function main()
	s = ArgParseSettings(
		description = "Plot isotropic Jij vs distance for all pairs. Optionally filter by element pair.",
	)
	@add_arg_table s begin
		"input"
		help = "Input XML file"
		required = true
		arg_type = String

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
	end
	args = parse_args(ARGS, s)

	# Validate: if one element is specified, both must be
	e1, e2 = args["element1"], args["element2"]
	if (e1 !== nothing && e2 === nothing) || (e1 === nothing && e2 !== nothing)
		error("--element1 and --element2 must be specified together")
	end

	distance_jij_pairs = collect_jij_all(args["input"]; element1 = e1, element2 = e2)

	plot_jij(
		distance_jij_pairs;
		invert_sign = args["invert-sign"],
		half_jij = args["half-jij"],
		ymin = args["ymin"],
		ymax = args["ymax"],
	)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
