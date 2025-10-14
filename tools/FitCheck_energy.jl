#!/usr/bin/env julia
using ArgParse
using Plots
using Printf
using Statistics

const MARKER_SIZE = 5	# marker size for scatter
const MARKER_ALPHA = 0.8	# marker transparency
const AXIS_PADDING = 10.0	# padding around data range (meV)

"""
	calculate_statistics(y1::Vector{Float64}, y2::Vector{Float64}) -> Dict{String, Float64}

	Compute summary statistics between observed and predicted values.
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
	parse_file(file::String) -> Tuple{Vector{Float64}, Vector{Float64}}

	Parse a text file and return observed and predicted vectors.
	- Ignore comment lines starting with `#` and empty lines
	- Allow either 2 columns (observed predicted) or >=3 columns (use 2nd and 3rd)
"""
function parse_file(file::String)::Tuple{Vector{Float64}, Vector{Float64}}
	raw_lines = readlines(file)
	data = Vector{NTuple{2, Float64}}()
	for (li, line) in enumerate(raw_lines)
		line_str = strip(line)
		isempty(line_str) && continue
		startswith(line_str, "#") && continue
		cols = split(line_str)
		if length(cols) == 2
			o = parse(Float64, cols[1])
			p = parse(Float64, cols[2])
			push!(data, (o, p))
		elseif length(cols) >= 3
			# For >=3 columns, use the 2nd and 3rd columns (consistent with yyplot)
			o = parse(Float64, cols[2])
			p = parse(Float64, cols[3])
			push!(data, (o, p))
		else
			error("Line $(li) has $(length(cols)) elements, expected 2 or >=3: $(line)")
		end
	end

	observed = Float64[o for (o, _) in data]
	predicted = Float64[p for (_, p) in data]
	return observed, predicted
end

"""
	create_plot() -> Tuple{Plots.Plot, String}

	Create a Plots.jl plot configured for energy comparison.
"""

function create_plot()::Tuple{Plots.Plot, String}
	title = "Energy Comparison"
	xlabel = "Observed Energy (meV)"
	ylabel = "Predicted Energy (meV)"
	unit = "meV"

	p = plot(title = title,
		xlabel = xlabel,
		ylabel = ylabel,
		legend = :topleft,
		grid = true,
		size = (1000, 1000),
		dpi = 300,
		framestyle = :box,
		aspect_ratio = 1)

	# Reference line y = x (wide enough, will be hidden by x/y limits later)
	plot!(p, [-100000.0, 100000.0], [-100000.0, 100000.0],
		line = (:black, 1), label = "y = x")

	return p, unit
end

"""
	plot_energy(files::Vector{String}; output::Union{String,Nothing}=nothing, lim::Union{Float64,Nothing}=nothing)

Load energy data from files (convert eV -> meV),
shift each file by its own observed-value center, then draw scatter and y=x line.
"""
function plot_energy(
	files::Vector{String};
	output::Union{String, Nothing} = nothing,
	lim::Union{Float64, Nothing} = nothing,
)
	observed_lists = Vector{Vector{Float64}}()
	predicted_lists = Vector{Vector{Float64}}()

	for file in files
		obs, pre = parse_file(file)
		obs .*= 1000.0	# eV -> meV
		pre .*= 1000.0
		push!(observed_lists, obs)
		push!(predicted_lists, pre)
	end

	# Shift origin per file by its own observed-value center
	for i in eachindex(observed_lists)
		center_i = (minimum(observed_lists[i]) + maximum(observed_lists[i])) / 2
		observed_lists[i] .-= center_i
		predicted_lists[i] .-= center_i
	end

	# Print statistics per dataset
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

	# Create base plot and add reference line
	p, unit = create_plot()

	# Determine axis limits
	if lim !== nothing
		xmin, xmax = -lim, lim
		ymin, ymax = -lim, lim
	else
		minv = min(minimum(vcat(observed_lists...)), minimum(vcat(predicted_lists...)))
		maxv = max(maximum(vcat(observed_lists...)), maximum(vcat(predicted_lists...)))
		xmin, xmax = minv - AXIS_PADDING, maxv + AXIS_PADDING
		ymin, ymax = minv - AXIS_PADDING, maxv + AXIS_PADDING
	end

	xlims!(p, (xmin, xmax))
	ylims!(p, (ymin, ymax))

	# Scatter for each file dataset
	for (i, (obs, pre)) in enumerate(zip(observed_lists, predicted_lists))
		series_label = @sprintf("%d: %s (RMSE: %.2f %s)", i, basename(files[i]), calculate_statistics(obs, pre)["RMSE"], unit)
		plot!(p, obs, pre,
			seriestype = :scatter,
			markersize = MARKER_SIZE,
			markeralpha = MARKER_ALPHA,
			label = series_label,
		)
	end

	if !isnothing(output)
		savefig(p, output)
		println("\nPlot saved as '$(output)'")
	else
		println("\nPlot not saved.")
	end

	display(p)
	println("Press Enter to exit...")
	readline()
	return nothing
end

# ---- CLI ----
s = ArgParseSettings(
	description = """
	Draw a meV scatter plot between observed and predicted energies (given in eV),\n
		after shifting each file by its own observed-value center. Input is whitespace-separated text.\n
		Accepts either 2 columns (observed predicted) or 3+ columns (index observed predicted).\n
		If --lim is provided, axes are fixed to [-lim, lim] (meV).
	""",
	version = "1.0.0",
	add_version = true,
)

@add_arg_table s begin
	"files"
	help = "Input data files (multiple files allowed)"
	arg_type = String
	nargs = '+'
	required = true

	"--output" "-o"
	help = "Output filename (format inferred from extension)"
	arg_type = String
	default = nothing

	"--lim" "-l"
	help = "X/Y axis limits (meV). Fix to [-lim, lim]"
	arg_type = Float64
	default = nothing
end

parsed_args = parse_args(s)
plot_energy(parsed_args["files"];
	output = parsed_args["output"],
	lim = parsed_args["lim"],
)


