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

# function to plot the histogram of magmom
function plot_histogram(
	input::AbstractString,
	n_atoms::Integer,
	target_atom_indices::AbstractVector{<:Integer},
	min_bound::Union{<:Real, Nothing},
	max_bound::Union{<:Real, Nothing},
	bin_width::Union{<:Real, Nothing},
	vertical_lines::Union{AbstractVector{<:Real}, Nothing},
)
	# check the input file
	if !isfile(input)
		error("The input file does not exist: $input")
	end

	# check the total atoms
	if n_atoms <= 0
		error("The total number of atoms must be positive: $n_atoms")
	end

	# read the input file
	embset::Vector{SpinConfig} = read_embset(input, n_atoms)

	# determine the target atoms
	target_atoms = (-1 in target_atom_indices) ? collect(1:n_atoms) : target_atom_indices

	# collect the magmom size
	magmom_size_list = Float64[
		sconfig.magmom_size[atom]
		for sconfig in embset
		for atom in target_atoms
	]

	# statistical analysis
	mean_magmom_size = mean(magmom_size_list)
	std_magmom_size = std(magmom_size_list)
	var_magmom_size = var(magmom_size_list)
	println("statistics of magnetic moment size:")
	println("mean: ", round(mean_magmom_size, digits = 4))
	println("variance: ", round(var_magmom_size, digits = 4))
	println("std: ", round(std_magmom_size, digits = 4))
	
	# set the bins parameter
	min_bound_value = something(min_bound, 0.0)
	max_bound_value = something(max_bound, maximum(magmom_size_list) + 1)
	bin_width_value = isnothing(bin_width) ? fd_bin_width(magmom_size_list) : float(bin_width)
	# ensure positive step; if not, fallback
	if !(bin_width_value > 0)
		bin_width_value = fd_bin_width(magmom_size_list)
	end
	println("bin width: ", round(bin_width_value, digits = 6))

	bins = range(min_bound_value, max_bound_value, step = bin_width_value)

	# plot the histogram
	p = histogram(
		magmom_size_list,
		bins = bins,
		xlabel = "Size of magnetic moment",
		ylabel = "Frequency",
		legend = false,
	)

	# plot the vertical lines
	if !isnothing(vertical_lines)
		for v::Float64 in vertical_lines
			vline!(p, [v], color = :black)
		end
	end

	display(p)
	println("\nPress Enter to exit the program.")
	readline()
end


s = ArgParseSettings(
	description = "Plot the histogram of magmom from an EMBSET.txt format file",
)
@add_arg_table s begin
	"--input", "-i"
	help = "The input file (i.e. EMBSET.txt)"
	required = true

	"--n_atoms", "-n"
	help = "The total number of atoms"
	required = true
	arg_type = Int

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
	parsed_args["input"],
	parsed_args["n_atoms"],
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
