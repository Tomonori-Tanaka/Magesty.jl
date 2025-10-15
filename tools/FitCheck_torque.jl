#!/usr/bin/env julia
using ArgParse
using Plots
using Printf
using Statistics
using LinearAlgebra

const MARKER_SIZE = 5# marker size for scatter
const MARKER_ALPHA = 0.8# marker transparency
const AXIS_PADDING = 10.0# padding around data range (meV)

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
	parse_torque_file(file::String) -> Vector{NamedTuple}

	Parse a torque file and return vector of torque data.
	Format: atom_index element DFT_torque_x DFT_torque_y DFT_torque_z SCE_torque_x SCE_torque_y SCE_torque_z
	- Ignore comment lines starting with `#` and empty lines
"""
function parse_torque_file(file::String)::Vector{NamedTuple}
	raw_lines = readlines(file)
	data = Vector{NamedTuple}()

	for (li, line) in enumerate(raw_lines)
		line_str = strip(line)
		isempty(line_str) && continue
		startswith(line_str, "#") && continue

		cols = split(line_str)
		if length(cols) >= 8
			atom_index = parse(Int, cols[1])
			element = cols[2]
			dft_x = parse(Float64, cols[3])
			dft_y = parse(Float64, cols[4])
			dft_z = parse(Float64, cols[5])
			sce_x = parse(Float64, cols[6])
			sce_y = parse(Float64, cols[7])
			sce_z = parse(Float64, cols[8])

			push!(
				data,
				(
					atom_index = atom_index,
					element = element,
					dft_torque = [dft_x, dft_y, dft_z],
					sce_torque = [sce_x, sce_y, sce_z],
				),
			)
		else
			error("Line $(li) has $(length(cols)) elements, expected >=8: $(line)")
		end
	end

	return data
end

"""
	filter_data(data::Vector{NamedTuple}, atom_indices::Union{Vector{Int}, Nothing}, elements::Union{Vector{String}, Nothing}) -> Vector{NamedTuple}

	Filter torque data by atom indices or elements.
	atom_indices and elements are mutually exclusive.
"""
function filter_data(
	data::Vector{NamedTuple},
	atom_indices::Union{Vector{Int}, Nothing},
	elements::Union{Vector{String}, Nothing},
)::Vector{NamedTuple}
	if atom_indices !== nothing && elements !== nothing
		error("atom_indices and elements options are mutually exclusive")
	end

	if atom_indices !== nothing
		return filter(x -> x.atom_index in atom_indices, data)
	elseif elements !== nothing
		return filter(x -> x.element in elements, data)
	else
		return data
	end
end

"""
	create_plot(plot_type::String) -> Tuple{Plots.Plot, String}

	Create a Plots.jl plot configured for torque comparison.
"""
function create_plot(plot_type::String)::Tuple{Plots.Plot, String}
	if plot_type == "all"
		title = "Torque Component Comparison"
		xlabel = "DFT Torque (meV)"
		ylabel = "SCE Torque (meV)"
		unit = "meV"
	elseif plot_type == "norm"
		title = "Torque Magnitude Comparison"
		xlabel = "DFT Torque Magnitude (meV)"
		ylabel = "SCE Torque Magnitude (meV)"
		unit = "meV"
	elseif plot_type == "dir"
		title = "Torque Direction Comparison"
		xlabel = "DFT Torque Direction Component"
		ylabel = "SCE Torque Direction Component"
		unit = "unitless"
	else
		error("Invalid plot_type: $plot_type")
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

	# Reference line y = x
	plot!(p, [-100000.0, 100000.0], [-100000.0, 100000.0],
		line = (:black, 1), label = "y = x")

	return p, unit
end

"""
	plot_torque(files::Vector{String}, plot_type::String; output::Union{String,Nothing}=nothing, lim::Union{Float64,Nothing}=nothing, atom_indices::Union{Vector{Int},Nothing}=nothing, elements::Union{Vector{String},Nothing}=nothing)

	Load torque data from files (convert eV -> meV),
	shift each file by its own observed-value center, then draw scatter and y=x line.
"""
function plot_torque(
	files::Vector{String},
	plot_type::String;
	output::Union{String, Nothing} = nothing,
	lim::Union{Float64, Nothing} = nothing,
	atom_indices::Union{Vector{Int}, Nothing} = nothing,
	elements::Union{Vector{String}, Nothing} = nothing,
)
	observed_lists = Vector{Vector{Float64}}()
	predicted_lists = Vector{Vector{Float64}}()

	for file in files
		data = parse_torque_file(file)
		filtered_data = filter_data(data, atom_indices, elements)

		if isempty(filtered_data)
			println("Warning: No data found in $file after filtering")
			continue
		end

		observed = Float64[]
		predicted = Float64[]

		if plot_type == "all"
			# Plot all x, y, z components
			for item in filtered_data
				obs = item.dft_torque .* 1000.0  # eV -> meV
				pre = item.sce_torque .* 1000.0
				append!(observed, obs)
				append!(predicted, pre)
			end
		elseif plot_type == "norm"
			# Plot magnitude
			for item in filtered_data
				obs_norm = norm(item.dft_torque) * 1000.0  # eV -> meV
				pre_norm = norm(item.sce_torque) * 1000.0
				push!(observed, obs_norm)
				push!(predicted, pre_norm)
			end
		elseif plot_type == "dir"
			# Plot direction components (unit vectors)
			for item in filtered_data
				obs_vec = item.dft_torque
				pre_vec = item.sce_torque
				obs_norm = norm(obs_vec)
				pre_norm = norm(pre_vec)

				if obs_norm > 1e-10 && pre_norm > 1e-10
					obs_unit = obs_vec ./ obs_norm
					pre_unit = pre_vec ./ pre_norm
					append!(observed, obs_unit)
					append!(predicted, pre_unit)
				end
			end
		end

		push!(observed_lists, observed)
		push!(predicted_lists, predicted)
	end

	if isempty(observed_lists)
		error("No data found in any file after filtering")
	end

	# Create base plot and add reference line
	p, unit = create_plot(plot_type)

	# Print statistics per dataset
	for (i, file) in enumerate(files)
		if !isempty(observed_lists[i])
			stats = calculate_statistics(observed_lists[i], predicted_lists[i])
			println("-"^30)
			println("File Index: $i")
			println("File Name: $file")
			println("RMSE: $(@sprintf("%.4f", stats["RMSE"])) $unit")
			println("R²: $(@sprintf("%.4f", stats["R²"]))")
			println("Max Error: $(@sprintf("%.4f", stats["Max Error"])) $unit")
			println("-"^30)
		end
	end

	# Determine axis limits
	if lim !== nothing
		xmin, xmax = -lim, lim
		ymin, ymax = -lim, lim
	else
		all_observed = vcat(observed_lists...)
		all_predicted = vcat(predicted_lists...)
		minv = min(minimum(all_observed), minimum(all_predicted))
		maxv = max(maximum(all_observed), maximum(all_predicted))
		xmin, xmax = minv - AXIS_PADDING, maxv + AXIS_PADDING
		ymin, ymax = minv - AXIS_PADDING, maxv + AXIS_PADDING
	end

	xlims!(p, (xmin, xmax))
	ylims!(p, (ymin, ymax))

	# Scatter for each file dataset
	for (i, (obs, pre)) in enumerate(zip(observed_lists, predicted_lists))
		if !isempty(obs)
			series_label = @sprintf(
				"%d: %s (RMSE: %.2f %s)",
				i,
				basename(files[i]),
				calculate_statistics(obs, pre)["RMSE"],
				unit
			)
			plot!(p, obs, pre,
				seriestype = :scatter,
				markersize = MARKER_SIZE,
				markeralpha = MARKER_ALPHA,
				label = series_label,
			)
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
		Draw a scatter plot between observed and predicted magnetic torque (given in eV),\n
		after shifting each file by its own observed-value center. Input is whitespace-separated text.\n
		Format: atom_index element DFT_torque_x DFT_torque_y DFT_torque_z SCE_torque_x SCE_torque_y SCE_torque_z\n
		If --lim is provided, axes are fixed to [-lim, lim] (meV).\n
		Use --atom-indices or --elements to filter data (mutually exclusive).
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

	"--type", "-t"
	help = "Type of data to plot (all: plot all torque components, norm: plot norm of torque, dir: plot direction of torque)"
	arg_type = String
	default = "all"
	range_tester = x -> x ∈ ["all", "norm", "dir"]

	"--output", "-o"
	help = "Output filename (png, svg, pdf; format inferred from extension)"
	arg_type = String
	default = nothing

	"--lim", "-l"
	help = "X/Y axis limits (meV). Fix to [-lim, lim]"
	arg_type = Float64
	default = nothing

	"--atom-indices", "-a"
	help = "Atom indices to plot (comma-separated, e.g., '1,2,3')"
	arg_type = String
	default = nothing

	"--elements", "-e"
	help = "Elements to plot (comma-separated, e.g., 'Fe,Co')"
	arg_type = String
	default = nothing
end

parsed_args = parse_args(s)

# Parse atom indices and elements
atom_indices = nothing
elements = nothing

if parsed_args["atom-indices"] !== nothing
	atom_indices = [parse(Int, x) for x in split(parsed_args["atom-indices"], ",")]
end

if parsed_args["elements"] !== nothing
	elements = String.(strip.(split(parsed_args["elements"], ",")))
end

plot_type = parsed_args["type"] === nothing ? "all" : parsed_args["type"]

println("Plotting $(plot_type) type of data")
println("Atom indices: $atom_indices")
println("Elements: $elements")
println("Output: $(parsed_args["output"])")
println("Lim: $(parsed_args["lim"])")


plot_torque(parsed_args["files"], plot_type;
	output = parsed_args["output"],
	lim = parsed_args["lim"],
	atom_indices = atom_indices,
	elements = elements,
)