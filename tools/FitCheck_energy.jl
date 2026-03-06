#!/usr/bin/env julia
using ArgParse
using Plots
using Printf
using Statistics

const MARKER_SIZE = 5	# marker size for scatter
const MARKER_ALPHA = 0.8	# marker transparency
const AXIS_PADDING = 10.0	# padding around data range (meV)

"""
	parse_index_list(spec::AbstractString) -> Vector{Int}

Parse a string like \"1,3,5-10\" into a list of integer indices.
Supports comma-separated integers and ranges with hyphens.
"""
function parse_index_list(spec::AbstractString)::Vector{Int}
	indices = Int[]
	for token in split(spec, ',')
		t = strip(token)
		if isempty(t)
			continue
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
	return indices
end

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
	xlabel = "DFT Energy (meV)"
	ylabel = "SCE Energy (meV)"
	unit = "meV"

	p = plot(title = title,
		xlabel = xlabel,
		ylabel = ylabel,
		legend = :topleft,
		grid = false,
		size = (1000, 1000),
		dpi = 300,
		framestyle = :box,
		aspect_ratio = 1)

	# Reference line y = x (wide enough, will be hidden by x/y limits later)
	plot!(p, [-100000.0, 100000.0], [-100000.0, 100000.0],
		line = (:black, 1), label = "y = x")
	# Dashed lines at x=0 and y=0 (GR: line = (width, :dash))
	plot!(p, [0.0, 0.0], [-100000.0, 100000.0];
		color = :gray, line = (1, :dash), label = "")
	plot!(p, [-100000.0, 100000.0], [0.0, 0.0];
		color = :gray, line = (1, :dash), label = "")

	return p, unit
end

"""
	plot_energy(files::Vector{String};
		output::Union{String,Nothing}=nothing,
		lim::Union{Float64,Nothing}=nothing,
		lim_min::Union{Float64,Nothing}=nothing,
		lim_max::Union{Float64,Nothing}=nothing,
		zero_min::Bool=false,
		zero_at::Union{Float64,Nothing}=nothing)

	Load energy data from files (convert eV -> meV),
	shift each file by its own observed-value center (or optionally so that the global minimum
	observed energy or a specified energy becomes 0), then draw scatter and y=x line.

	If `lim` is provided, X/Y axes are fixed to [-lim, lim] (meV).
	If `lim_min` and/or `lim_max` are provided, X/Y axes are fixed to [lim_min, lim_max] (meV).
	When only one of `lim_min` or `lim_max` is given, the other is set to the same absolute value
	with the opposite sign.
	If `zero_min` is true, all datasets are shifted so that the global minimum observed energy
	becomes 0.
	If `zero_at` is provided (in meV), all datasets are shifted so that the specified energy
	corresponds to 0. `zero_min` and `zero_at` are mutually exclusive.
"""
function plot_energy(
	files::Vector{String};
	output::Union{String, Nothing} = nothing,
	lim::Union{Float64, Nothing} = nothing,
	lim_min::Union{Float64, Nothing} = nothing,
	lim_max::Union{Float64, Nothing} = nothing,
	zero_min::Bool = false,
	zero_at::Union{Float64, Nothing} = nothing,
	colored_indices::Union{Vector{Int}, Nothing} = nothing,
	no_legend::Bool = false,
	marker_size::Real = MARKER_SIZE,
	marker_alpha::Real = MARKER_ALPHA,
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

	if zero_min && zero_at !== nothing
		error("Options `--zero-min` and `--zero-at` are mutually exclusive.")
	end

	if zero_at !== nothing
		# Shift all datasets so that the specified energy becomes 0
		for i in eachindex(observed_lists)
			observed_lists[i] .-= zero_at
			predicted_lists[i] .-= zero_at
		end
	elseif zero_min
		# Shift all datasets so that the global minimum observed energy becomes 0
		global_min = minimum(vcat(observed_lists...))
		for i in eachindex(observed_lists)
			observed_lists[i] .-= global_min
			predicted_lists[i] .-= global_min
		end
	else
		# Shift origin per file by its own observed-value center
		for i in eachindex(observed_lists)
			center_i = (minimum(observed_lists[i]) + maximum(observed_lists[i])) / 2
			observed_lists[i] .-= center_i
			predicted_lists[i] .-= center_i
		end
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

	# Optionally hide legend
	if no_legend
		plot!(p, legend = false)
	end

	# Determine axis limits
	if lim !== nothing
		xmin, xmax = -lim, lim
		ymin, ymax = -lim, lim
	elseif lim_min !== nothing || lim_max !== nothing
		local_min = lim_min
		local_max = lim_max
		if local_min === nothing && local_max !== nothing
			local_min = -abs(local_max)
		elseif local_min !== nothing && local_max === nothing
			local_max = abs(local_min)
		end
		# At this point both local_min and local_max must be non-nothing
		xmin, xmax = local_min, local_max
		ymin, ymax = local_min, local_max
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
		series_label = @sprintf(
			"%d: %s (RMSE: %.2f %s)",
			i,
			basename(files[i]),
			calculate_statistics(obs, pre)["RMSE"],
			unit,
		)

		if colored_indices === nothing
			# Plot all points with the same marker style
			plot!(p, obs, pre,
				seriestype = :scatter,
				markersize = marker_size,
				markeralpha = marker_alpha,
				label = series_label,
			)
		else
			n = length(obs)
			mask_colored = falses(n)
			for idx in colored_indices
				if 1 <= idx <= n
					mask_colored[idx] = true
				end
			end
			mask_uncolored = .!mask_colored

			obs_uncolored = obs[mask_uncolored]
			pre_uncolored = pre[mask_uncolored]
			obs_colored = obs[mask_colored]
			pre_colored = pre[mask_colored]

			if !isempty(obs_uncolored)
				plot!(p, obs_uncolored, pre_uncolored,
					seriestype = :scatter,
					markersize = marker_size,
					markeralpha = marker_alpha,
					label = series_label,
				)
			end
			if !isempty(obs_colored)
				plot!(p, obs_colored, pre_colored,
					seriestype = :scatter,
					markersize = marker_size,
					markeralpha = marker_alpha,
					color = :deepskyblue,
					label = string(series_label, " (colored)"),
				)
			end
		end
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
		after shifting each file by its own observed-value center (or optionally so that the global minimum observed energy becomes 0, or a specified energy becomes 0). Input is whitespace-separated text.\n
		Accepts either 2 columns (observed predicted) or 3+ columns (index observed predicted).\n
		If --lim is provided, axes are fixed to [-lim, lim] where lim is specified in eV (internally converted to meV).\n
		Alternatively, use --lim-min/--lim-max to fix axes to [lim-min, lim-max], also specified in eV.
	""",
	version = "0.1.0",
	add_version = true,
)

@add_arg_table s begin
	"files"
	help = "Input data files (multiple files allowed)"
	arg_type = String
	nargs = '+'
	required = true

	"--output", "-o"
	help = "Output filename (png, svg, pdf; format inferred from extension)"
	arg_type = String
	default = nothing

	"--lim", "-l"
	help = "X/Y axis limits (eV). Fix to [-lim, lim] (converted to meV internally)"
	arg_type = Float64
	default = nothing

	"--lim-min"
	help = "Minimum X/Y axis limit (eV). Used together with --lim-max; if one side is omitted, it is set to the opposite sign with the same absolute value."
	arg_type = Float64
	default = nothing

	"--lim-max"
	help = "Maximum X/Y axis limit (eV). Used together with --lim-min; if one side is omitted, it is set to the opposite sign with the same absolute value."
	arg_type = Float64
	default = nothing

	"--zero-min", "-z"
	help = "Shift all energies so that the global minimum observed energy becomes 0"
	action = :store_true

	"--zero-at"
	help = "Shift all energies so that the specified observed energy (eV) becomes 0 (mutually exclusive with --zero-min; converted to meV internally)"
	arg_type = Float64
	default = nothing

	"--colored", "-c"
	help = "Indices to highlight in a different color (e.g., '5,7-10')"
	arg_type = String
	default = nothing

	"--marker-size", "-m"
	help = "Marker size for scatter points"
	arg_type = Float64
	default = MARKER_SIZE

	"--marker-alpha", "-A"
	help = "Marker transparency (0-1) for scatter points"
	arg_type = Float64
	default = MARKER_ALPHA

	"--no-legend", "-L"
	help = "Disable legend in the plot"
	action = :store_true
end

parsed_args = parse_args(s)

# Convert CLI energy options from eV to meV for internal use
lim_mev = parsed_args["lim"]
if lim_mev !== nothing
	lim_mev *= 1000.0
end

lim_min_mev = parsed_args["lim-min"]
if lim_min_mev !== nothing
	lim_min_mev *= 1000.0
end

lim_max_mev = parsed_args["lim-max"]
if lim_max_mev !== nothing
	lim_max_mev *= 1000.0
end

zero_at_mev = parsed_args["zero-at"]
if zero_at_mev !== nothing
	zero_at_mev *= 1000.0
end

colored_indices = nothing
if parsed_args["colored"] !== nothing
	colored_indices = parse_index_list(parsed_args["colored"])
end

plot_energy(parsed_args["files"];
	output = parsed_args["output"],
	lim = lim_mev,
	lim_min = lim_min_mev,
	lim_max = lim_max_mev,
	zero_min = parsed_args["zero-min"],
	zero_at = zero_at_mev,
	colored_indices = colored_indices,
	no_legend = parsed_args["no-legend"],
	marker_size = parsed_args["marker-size"],
	marker_alpha = parsed_args["marker-alpha"],
)


