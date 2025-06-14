#!/usr/bin/env julia
using ArgParse
using Plots
using Printf
using LinearAlgebra
using Statistics

"""
	create_plot(type::String) -> Plots.Plot

Create a plot with appropriate settings based on the data type.

# Arguments
- `type::String`: Type of data ('e' for energy, 'm' for magnetic field)

# Returns
A Plots.Plot object with appropriate settings
"""
function create_plot(type::String)::Tuple{Plots.Plot, String}
	if type == "e"
		title = "Energy Comparison"
		xlabel = "Observed Energy (meV)"
		ylabel = "Predicted Energy (meV)"
		unit = "meV"
	elseif type == "m"
		title = "Magnetic Field Comparison"
		xlabel = "Observed Magnetic Field (meV)"
		ylabel = "Predicted Magnetic Field (meV)"
		unit = "meV"
	else
		error("Invalid type: $type. Must be 'e' or 'm'")
	end

	p = plot(title = title,
		xlabel = xlabel,
		ylabel = ylabel,
		legend = :topleft,
		grid = true,
		size = (1000, 1000),
		dpi = 300,
		framestyle = :box,
		aspect_ratio = 1)

	# Plot the reference line (y = x)
	plot!(p, [-100000, 100000], [-100000, 100000],
		line = (:black, 1),
		label = "y = x")

	return p, unit
end

"""
	calculate_statistics(y1::Vector{Float64}, y2::Vector{Float64}) -> Dict{String, Float64}

Calculate statistics for the given data.

# Arguments
- `y1::Vector{Float64}`: Observed values
- `y2::Vector{Float64}`: Predicted values

# Returns
Dictionary containing various statistics
"""
function calculate_statistics(y1::Vector{Float64}, y2::Vector{Float64})::Dict{String, Float64}
	rmse = sqrt(mean((y1 .- y2) .^ 2))
	max_error = maximum(abs.(y1 .- y2))
	mean_error = mean(abs.(y1 .- y2))
	std_error = std(y1 .- y2)
	r_squared = 1 - sum((y1 .- y2) .^ 2) / sum((y1 .- mean(y1)) .^ 2)

	return Dict(
		"RMSE" => rmse,
		"Max Error" => max_error,
		"Mean Error" => mean_error,
		"Std Error" => std_error,
		"R²" => r_squared,
	)
end

"""
	parse_file(file::String)::Tuple{Vector{Float64}, Vector{Float64}}

Parse the file and return the observed and predicted values.
Comment lines (starting with #, regardless of the position) and empty lines are skipped.

Assumed file format:
```
# Observed_Energy, Predicted_Energy
 1.0000000000, 1.0000000000
 2.0000000000, 2.0000000000
 3.0000000000, 3.0000000000
```

# Arguments
- `file::String`: Input file path

"""
function parse_file(file::String)::Tuple{Vector{Float64}, Vector{Float64}}
	# Read the file and filter out comment lines
	raw_lines = readlines(file)

	# Parse the data. Skip comment lines and empty lines.
	parsed_lines::Vector{Vector{Float64}} =
		[
			parse.(Float64, split(line)) for
			line in raw_lines if !isempty(strip(line)) && !startswith(strip(line), "#")
		]

	# Validate that each line has exactly two elements
	for (i, line) in enumerate(parsed_lines)
		if length(line) != 2
			error("Line $i has $(length(line)) elements, expected 2")
		end
	end

	# Convert to two vectors
	observed = [line[1] for line in parsed_lines]
	predicted = [line[2] for line in parsed_lines]

	return observed, predicted
end

"""
	main(files::Vector{String}; output::Union{String,Nothing}=nothing, lim::Union{Float64,Nothing}=nothing)

Plots energy data from input files by comparing observed and predicted values.

# Arguments
- `files::Vector{String}`: List of input file paths containing energy data.
- `output::Union{String,Nothing}`: Output file name for the plot. The PNG file is saved only if this argument is provided.
- `lim::Union{Float64,Nothing}`: Axis limits for both X and Y axes. If not provided, the limits will be computed from the data.
"""
function main(
	files::Vector{String},
	type::String,
	;
	output::Union{String, Nothing} = nothing,
	lim::Union{Float64, Nothing} = nothing)

	# result lists of list to store observed and predicted values, and those will be used for plotting
	observed_lists = Vector{Vector{Float64}}()
	predicted_lists = Vector{Vector{Float64}}()

	for file in files
		# Parse the file and return the observed and predicted values
		observed, predicted = parse_file(file)
		# Convert to meV from eV
		observed *= 1000
		predicted *= 1000

		push!(observed_lists, observed)
		push!(predicted_lists, predicted)
	end

	# Shift the all values by the center of the observed values if type is "e" (energy)
	if type == "e"
		center_observed = (minimum(vcat(observed_lists...)) + maximum(vcat(observed_lists...))) / 2
		for list in observed_lists
			list .-= center_observed
		end
		for list in predicted_lists
			list .-= center_observed
		end
	end

	# Print the statistics in each file data set 
	for (i, file) in enumerate(files)
		stats = calculate_statistics(observed_lists[i], predicted_lists[i])
		println("-"^30)
		println("File Index: $i")
		println("File Name: $file")
		println("RMSE: $(@sprintf("%.4f", stats["RMSE"])) meV")
		println("R²: $(@sprintf("%.4f", stats["R²"]))")
		println("Max Error: $(@sprintf("%.4f", stats["Max Error"])) meV")
		println("-"^30)
	end

	# Create the plot
	p, unit = create_plot(type)
	for (i, (observed, predicted)) in enumerate(zip(observed_lists, predicted_lists))
		plot!(p, observed, predicted,
			seriestype = :scatter,
			markersize = 5,
			markeralpha = 0.5,
			label = "$i: $(basename(files[i])) (RMSE: $(@sprintf("%.2f", calculate_statistics(observed, predicted)["RMSE"])) $unit)",
		)
	end

	# Set axis limits based on user input or overall data
	if lim !== nothing
		xlim_val = (-lim, lim)
		ylim_val = (-lim, lim)
	else
		minimum_value = min(minimum(vcat(observed_lists...)), minimum(vcat(predicted_lists...)))
		maximum_value = max(maximum(vcat(observed_lists...)), maximum(vcat(predicted_lists...)))
		xlim_val = (minimum_value - 10, maximum_value + 10)
		ylim_val = (minimum_value - 10, maximum_value + 10)
	end
	xlims!(p, xlim_val)
	ylims!(p, ylim_val)

	# Save the plot only if the output argument is provided
	if !isnothing(output)
		savefig(p, output)
		println("\nPlot saved as '$output'")
	else
		println("\nPlot not saved.")
	end
	display(p)

	println("Press Enter to exit...")
	readline()
end


# Set up command-line argument parsing
s = ArgParseSettings(
	description = """
		Generates a comparison plot between observed and predicted energy or magnetic field values.\n
		Input files should be space-separated text files with at least two columns:\n
		Column 1: Observed value, Column 2: Predicted value.\n
		\n
		Note that the data value is shifted by the center of the observed values if type is "e"
		(the value of shift is common for all file data sets)\n
	""",
	version = "1.0.0",
	add_version = true,
)

@add_arg_table s begin
	"files"
	help = "Input data files. Multiple files are allowed."
	arg_type = String
	nargs = '+'
	required = true

	"--type", "-t"
	help = "Type of data to plot (e: energy, m: magnetic field)"
	arg_type = String
	action = :store_arg
	required = true
	range_tester = x -> x ∈ ["e", "m"]

	"--output", "-o"
	help = "Output file name for the plot (if specified, the PNG file will be saved)"
	arg_type = String
	default = nothing

	"--lim", "-l"
	help = "Axis limits for X and Y axes (min max)"
	arg_type = Float64
	default = nothing
end

# Parse arguments and convert the lim argument from Vector to Tuple if provided
parsed_args = parse_args(s)

main(parsed_args["files"],
	parsed_args["type"],
	;
	output = parsed_args["output"],
	lim = parsed_args["lim"],
)
