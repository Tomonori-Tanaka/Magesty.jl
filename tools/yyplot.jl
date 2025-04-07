#!/usr/bin/env julia
using ArgParse
using Plots
using Printf
using LinearAlgebra
using Statistics

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

	# Validate input files
	for file in files
		if !isfile(file)
			error("Error: File not found: $file")
		end
	end

	# Arrays to store combined data from all files
	global_y1 = Float64[]
	global_y2 = Float64[]

	# Create the plot
	p, unit = create_plot(type)

	# Dictionary to store statistics for each file
	stats = Dict{String, Dict{String, Float64}}()

	# Process each file
	for file in files
		try
			# Read the file and filter out comment lines
			raw_lines = readlines(file)
			parsed_lines =
				[parse.(Float64, split(line)) for line in raw_lines if !startswith(line, "#")]

			# Check if there is valid data
			if isempty(parsed_lines)
				error("No valid data found in file $file")
			end

			# Convert each row to a 1xN matrix and combine vertically
			data = reduce(vcat, [reshape(row, 1, :) for row in parsed_lines])

			# Ensure the data has at least three columns
			if size(data, 2) < 3
				error("File $file does not contain enough columns")
			end

			# Preprocess data according to type
			y1, y2 = preprocess_data(data, type)

			# Append data to global arrays for computing overall limits
			append!(global_y1, y1)
			append!(global_y2, y2)

			# Calculate statistics
			stats[basename(file)] = calculate_statistics(y1, y2)

			# Add the data as a scatter plot for the current file
			plot!(p, y1, y2,
				seriestype = :scatter,
				markersize = 6,
				markeralpha = 0.7,
				label = "$(basename(file)) (RMSE: $(@sprintf("%.2f", stats[basename(file)]["RMSE"])) $unit)")
		catch e
			@error "Error processing file $file" exception = (e, catch_backtrace())
			continue
		end
	end

	# Set axis limits based on user input or overall data
	if lim !== nothing
		xlim_val = (-lim, lim)
		ylim_val = (-lim, lim)
	else
		if isempty(global_y1)
			error("Error: No valid data could be loaded from any file")
		end
		xlim_val = (minimum(global_y1), maximum(global_y1))
		ylim_val = (minimum(global_y2), maximum(global_y2))
	end
	xlims!(p, xlim_val)
	ylims!(p, ylim_val)

	# Print statistics for each file
	println("\nStatistics:")
	println("="^50)
	for (file, stat) in stats
		println("\nFile: $file")
		println("-"^30)
		for (key, value) in stat
			println("  $key: $(@sprintf("%.2f", value)) $unit")
		end
	end
	println("="^50)

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

"""
	preprocess_data(data::Matrix{Float64}, type::String) -> Tuple{Vector{Float64}, Vector{Float64}}

Preprocess the data according to the specified type.

# Arguments
- `data::Matrix{Float64}`: Input data matrix
- `type::String`: Type of data ('e' for energy, 'm' for magnetic field)

# Returns
Tuple of (y1, y2) where:
- y1: Processed observed values
- y2: Processed predicted values
"""
function preprocess_data(data::Matrix{Float64}, type::String)
	if type == "e"
		# Extract observed energy (column 2) and predicted energy (column 3)
		y1 = data[:, 2] * 1000  # Convert to meV
		y2 = data[:, 3] * 1000  # Convert to meV

		# Shift values by the mean of y1
		y1_mean = mean(y1)
		y1 .-= y1_mean
		y2 .-= y1_mean
	elseif type == "m"
		# Extract observed and predicted magnetic field values
		y1 = data[:, 2]  # No unit conversion needed
		y2 = data[:, 3]  # No unit conversion needed
	else
		error("Invalid type: $type. Must be 'e' or 'm'")
	end

	return y1, y2
end

"""
	create_plot(type::String) -> Plots.Plot

Create a plot with appropriate settings based on the data type.

# Arguments
- `type::String`: Type of data ('e' for energy, 'm' for magnetic field)

# Returns
A Plots.Plot object with appropriate settings
"""
function create_plot(type::String)
	if type == "e"
		title = "Energy Comparison"
		xlabel = "Observed Energy (meV)"
		ylabel = "Predicted Energy (meV)"
		unit = "meV"
	elseif type == "m"
		title = "Magnetic Field Comparison"
		xlabel = "Observed Magnetic Field (T)"
		ylabel = "Predicted Magnetic Field (T)"
		unit = "T"
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
function calculate_statistics(y1::Vector{Float64}, y2::Vector{Float64})
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

# Set up command-line argument parsing
s = ArgParseSettings(
	description = """
		Generates a comparison plot between observed and predicted energy or magnetic field values.
		Input files should be space-separated text files with at least three columns:
		Column 1: Index, Column 2: Observed value, Column 3: Predicted value.
	""",
	version = "1.0.0",
	add_version = true,
)

@add_arg_table s begin
	"files"
	help = "Input data files (space-separated text files)"
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
