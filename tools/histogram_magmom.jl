using ArgParse
using Plots
using Statistics

include("../src/SpinConfig.jl")
using .SpinConfigs

# function to plot the histogram of magmom
function plot_histogram(
	input_path::AbstractString,
	n_atoms::Integer,
	target_atom_indices::AbstractVector{<:Integer},
	min_bound::Union{<:Real, Nothing},
	max_bound::Union{<:Real, Nothing},
	bin_width::Union{<:Real, Nothing},
	vertical_lines::Union{AbstractVector{<:Real}, Nothing},
)
	# check the input file
	if !isfile(input_path)
		error("The input file does not exist: $input_path")
	end

	# check the total atoms
	if n_atoms <= 0
		error("The total number of atoms must be positive: $n_atoms")
	end

	# read the input file
	embset::Vector{SpinConfig} = read_embset(input_path, n_atoms)

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
	bin_width_value = something(bin_width, 0.1)

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
	"--input_path", "-i"
	help = "The input file (i.e. EMBSET.txt)"
	required = true

	"--n_atoms", "-n"
	help = "The total number of atoms"
	required = true
	arg_type = Int

	"--atoms", "-a"
	help = "The atoms to plot. If -1 is given (default), all atoms are plotted."
	nargs = '+'
	default = [-1]
	arg_type = Int

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
	parsed_args["input_path"],
	parsed_args["n_atoms"],
	parsed_args["atoms"],
	parsed_args["min_bound"],
	parsed_args["max_bound"],
	parsed_args["bin_width"],
	parsed_args["vertical_lines"],
)
