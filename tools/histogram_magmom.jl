using ArgParse
using Plots
using Statistics

include("../src/SpinConfigs.jl")
using .SpinConfigs

# parse atoms argument that may include ranges like "1-5" and comma-separated tokens
function parse_atom_indices(atoms_args::AbstractVector{<:AbstractString})::Vector{Int}
	indices = Int[]
	for raw_token in atoms_args
		for token in split(raw_token, ',')
			t = strip(token)
			if isempty(t)
				continue
			end
			if t == "-1"
				return [-1]
			end
			if occursin('-', t)
				parts = split(t, '-')
				if length(parts) != 2
					error("Invalid range token: $t")
				end
				start_str, end_str = strip.(parts)
				start_idx = tryparse(Int, start_str)
				end_idx = tryparse(Int, end_str)
				if isnothing(start_idx) || isnothing(end_idx)
					error("Invalid integer in range: $t")
				end
				if start_idx > end_idx
					error("Range start greater than end: $t")
				end
				append!(indices, collect(start_idx:end_idx))
			else
				val = tryparse(Int, t)
				if isnothing(val)
					error("Invalid integer token: $t")
				end
				push!(indices, val)
			end
		end
	end
	return indices
end

# Freedman–Diaconis rule bin width with sensible fallbacks
function fd_bin_width(values::AbstractVector{<:Real})::Float64
	v = collect(skipmissing(values))
	if isempty(v)
		return 0.1
	end
	n = length(v)
	if n < 2
		return 0.1
	end
	q75 = quantile(v, 0.75)
	q25 = quantile(v, 0.25)
	iqr = q75 - q25
	if iqr > 0
		h = 2 * iqr / (n^(1/3))
		return h > 0 ? float(h) : 0.1
	end
	# fallback to Scott's rule if IQR == 0
	s = std(v)
	if s > 0
		h = 3.49 * s / (n^(1/3))
		return h > 0 ? float(h) : 0.1
	end
	# final fallback: split range into 10 bins
	rng = maximum(v) - minimum(v)
	return rng > 0 ? float(rng / 10) : 0.1
end

# function to plot the histogram of magmom for multiple files
function plot_histogram(
	inputs::AbstractVector{<:AbstractString},
	target_atom_indices::AbstractVector{<:Integer},
	min_bound::Union{<:Real, Nothing},
	max_bound::Union{<:Real, Nothing},
	bin_width::Union{<:Real, Nothing},
	vertical_lines::Union{AbstractVector{<:Real}, Nothing},
)
	# check the input files
	for input in inputs
		if !isfile(input)
			error("The input file does not exist: $input")
		end
	end

	# Detect number of atoms from the first file
	n_atoms = detect_num_atoms(inputs[1])
	if n_atoms <= 0
		error("The total number of atoms must be positive: $n_atoms")
	end

	# determine the target atoms
	target_atoms = (-1 in target_atom_indices) ? collect(1:n_atoms) : target_atom_indices
	# validate and normalize atom indices when not selecting all
	if !(-1 in target_atom_indices)
		for idx in target_atoms
			if idx < 1 || idx > n_atoms
				error("Atom index out of range (1..$n_atoms): $idx")
			end
		end
		# remove duplicates and sort for stable processing
		target_atoms = sort(unique(target_atoms))
	end

	# collect data from all files
	all_magmom_data = Vector{Float64}[]
	file_stats = []
	
	for (file_idx, input) in enumerate(inputs)
		# read the input file
		embset::Vector{SpinConfig} = SpinConfigs.read_embset(input)
		
		# collect the magmom size for this file
		magmom_size_list = Float64[
			sconfig.magmom_size[atom]
			for sconfig in embset
			for atom in target_atoms
		]
		
		push!(all_magmom_data, magmom_size_list)
		
		# statistical analysis for this file
		mean_magmom_size = mean(magmom_size_list)
		std_magmom_size = std(magmom_size_list)
		var_magmom_size = var(magmom_size_list)
		
		file_stat = (
			filename = basename(input),
			mean = mean_magmom_size,
			std = std_magmom_size,
			var = var_magmom_size,
			count = length(magmom_size_list)
		)
		push!(file_stats, file_stat)
		
		println("File $file_idx: $(file_stat.filename)")
		println("  statistics of magnetic moment size:")
		println("  mean: ", round(mean_magmom_size, digits = 8))
		println("  variance: ", round(var_magmom_size, digits = 8))
		println("  std: ", round(std_magmom_size, digits = 8))
		println("  count: ", file_stat.count)
	end

	# combine all data for overall statistics and bin calculation
	combined_magmom_data = vcat(all_magmom_data...)
	
	# overall statistical analysis
	overall_mean = mean(combined_magmom_data)
	overall_std = std(combined_magmom_data)
	overall_var = var(combined_magmom_data)
	println("\nOverall statistics (combined data):")
	println("mean: ", round(overall_mean, digits = 8))
	println("variance: ", round(overall_var, digits = 8))
	println("std: ", round(overall_std, digits = 8))
	println("total count: ", length(combined_magmom_data))
	
	# set the bins parameter using combined data
	min_bound_value = something(min_bound, 0.0)
	max_bound_value = something(max_bound, maximum(combined_magmom_data) + 1)
	bin_width_value = isnothing(bin_width) ? fd_bin_width(combined_magmom_data) : float(bin_width)
	# ensure positive step; if not, fallback
	if !(bin_width_value > 0)
		bin_width_value = fd_bin_width(combined_magmom_data)
	end
	# ensure valid range
	if max_bound_value <= min_bound_value
		max_bound_value = min_bound_value + bin_width_value * 10
	end
	println("bin width: ", round(bin_width_value, digits = 6))

	bins = range(min_bound_value, max_bound_value, step = bin_width_value)
	
	# compute mode using binned data for combined data
	bin_centers = collect(bins)
	bin_counts = zeros(Int, length(bin_centers))
	for val in combined_magmom_data
		# find the closest bin center
		bin_idx = argmin(abs.(bin_centers .- val))
		bin_counts[bin_idx] += 1
	end
	max_count = maximum(bin_counts)
	if max_count > 0
		mode_bins = bin_centers[bin_counts .== max_count]
		println("mode (binned): ", join(round.(mode_bins, digits = 8), ", "))
	else
		println("mode (binned): (none)")
	end

	# plot the histogram with multiple files
	p = plot()
	
	# Define colors for different files
	colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :olive, :cyan]
	
	for (file_idx, (magmom_data, file_stat)) in enumerate(zip(all_magmom_data, file_stats))
		color = colors[mod1(file_idx, length(colors))]
		histogram!(
			p,
			magmom_data,
			bins = bins,
			alpha = 0.6,
			label = "$(file_stat.filename) (n=$(file_stat.count))",
			color = color,
			linecolor = :black,
			linewidth = 0.5
		)
	end
	
	# Combined histogram is not plotted
	
	# Set labels and title
	xlabel!(p, "Size of magnetic moment")
	ylabel!(p, "Frequency")
	title!(p, "Magnetic Moment Magnitude Distribution")
	
	# plot the vertical lines
	if !isnothing(vertical_lines)
		for v in vertical_lines
			vline!(p, [v], color = :black, linestyle = :dash, linewidth = 2)
		end
	end

	display(p)
	println("\nPress Enter to exit the program.")
	readline()
end


s = ArgParseSettings(
	description = "Plot the histogram of magmom from EMBSET.txt format files (supports multiple files)",
)
@add_arg_table s begin
	"inputs"
	help = "The input files (i.e. EMBSET.txt). Multiple files can be specified."
	nargs = '+'
	required = true

	"--atoms", "-a"
	help = "The atoms to plot. Accepts integers, ranges like 1-5, and comma-separated lists. If -1 is given (default), all atoms are plotted."
	nargs = '+'
	default = ["-1"]
	arg_type = String

	"--min_bound", "-l"
	help = "The lower bound of the histogram"
	arg_type = Float64

	"--max_bound", "-u"
	help = "The upper bound of the histogram"
	arg_type = Float64

	"--bin_width", "-w"
	help = "The bin_width size of the histogram"
	arg_type = Float64

	"--vertical_lines", "-v"
	help = "The values of the vertical lines"
	nargs = '+'
	arg_type = Float64
end

parsed_args = parse_args(ARGS, s)

plot_histogram(
	String.(parsed_args["inputs"]),
	begin
		atoms_arg = parsed_args["atoms"]
		# atoms_arg is Vector{String}; convert to Vector{Int} expanding ranges
		parse_atom_indices(atoms_arg)
	end,
	parsed_args["min_bound"],
	parsed_args["max_bound"],
	parsed_args["bin_width"],
	parsed_args["vertical_lines"],
)
