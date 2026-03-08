#!/usr/bin/env julia
#
# MFA sampling with quantization axes along Lebedev grid directions.
# Grid: either generate in script with --order N (Lebedev.jl) or load from file with --lebedev FILE.
# Include sampling_mfa.jl to reuse MFA logic and VASP tools.
#

include("sampling_mfa.jl")
using Lebedev

"""
	load_lebedev_axes(path::AbstractString) -> (Vector{Vector{Float64}}, Union{Nothing,Vector{Float64}})

Load unit vectors (quantization axes) from a Lebedev grid file.
Format: one line per point, "x y z" or "x y z weight".
Returns (axes, weights); weights is nothing if no weight column.
"""
function load_lebedev_axes(path::AbstractString)
	axes = Vector{Vector{Float64}}()
	weights = Float64[]
	has_weight = false
	open(path) do f
		for line in eachline(f)
			line = strip(line)
			isempty(line) && continue
			tok = split(line)
			length(tok) < 3 && continue
			x, y, z = parse(Float64, tok[1]), parse(Float64, tok[2]), parse(Float64, tok[3])
			v = [x, y, z]
			n = norm(v)
			iszero(n) && continue
			push!(axes, v ./ n)
			if length(tok) >= 4
				push!(weights, parse(Float64, tok[4]))
				has_weight = true
			end
		end
	end
	isempty(axes) && error("No axes found in Lebedev file: $path")
	w = has_weight && length(weights) == length(axes) ? weights : nothing
	return axes, w
end

"""
	generate_lebedev_grid(order::Integer) -> (Vector{Vector{Float64}}, Vector{Float64})

Generate Lebedev grid in script using Lebedev.jl (order 3,5,7,...,131; e.g. order 9 -> 38 points).
Returns (axes, weights) on unit sphere.
"""
function generate_lebedev_grid(order::Integer)
	x, y, z, w = lebedev_by_order(order)
	axes = [[x[i], y[i], z[i]] for i in eachindex(x)]
	return axes, collect(w)
end

"""
	filter_axes_by_angle(axes, weights, theta_min_deg, theta_max_deg, phi_min_deg, phi_max_deg)

Filter axes (unit vectors x,y,z) by polar angle theta [0,180] and azimuth phi (deg).
theta = 0 is +z; phi in degrees, same convention as lebedev.py.
weights can be nothing; then returned weights are nothing.
"""
function filter_axes_by_angle(
	axes::Vector{Vector{Float64}},
	weights::Union{Nothing,Vector{Float64}},
	theta_min_deg::Float64,
	theta_max_deg::Float64,
	phi_min_deg::Float64,
	phi_max_deg::Float64,
)
	theta_min = deg2rad(theta_min_deg)
	theta_max = deg2rad(theta_max_deg)
	phi_min = deg2rad(phi_min_deg)
	phi_max = deg2rad(phi_max_deg)

	out_axes = Vector{Vector{Float64}}()
	out_weights = weights === nothing ? nothing : Float64[]

	for (k, v) in enumerate(axes)
		x, y, z = v[1], v[2], v[3]
		z = clamp(z, -1.0, 1.0)
		theta = acos(z)
		phi = atan(y, x)

		ok_theta = theta >= theta_min && theta <= theta_max
		ok_phi = true
		if phi_max - phi_min < 2 * π - 1e-9
			if phi_max >= phi_min
				ok_phi = phi >= phi_min && phi <= phi_max
			else
				ok_phi = (phi >= phi_min) || (phi <= phi_max)
			end
		end

		if ok_theta && ok_phi
			push!(out_axes, v)
			if out_weights !== nothing
				push!(out_weights, weights[k])
			end
		end
	end

	return out_axes, out_weights
end

"""
	write_grid_file(path, axes, weights)

Write grid points to file: one line per point "x y z" or "x y z weight".
"""
function write_grid_file(
	path::AbstractString,
	axes::Vector{Vector{Float64}},
	weights::Union{Nothing,Vector{Float64}},
)
	open(path, "w") do f
		for (i, v) in enumerate(axes)
			if weights !== nothing && i <= length(weights)
				println(f, v[1], " ", v[2], " ", v[3], " ", weights[i])
			else
				println(f, v[1], " ", v[2], " ", v[3])
			end
		end
	end
end

"""
	sampling_mfa_lebedev(input_path, value, axes, n_per_axis, output_start, digits, variable, fix_spec)

Generate MFA samples with quantization axes from the Lebedev grid.
n_per_axis: number of samples per Lebedev direction (total = length(axes) * n_per_axis).
output_start: base index for sample IDs (0 for a single run).
"""
function sampling_mfa_lebedev(
	input_path::AbstractString,
	value::Real,
	axes::Vector{Vector{Float64}},
	n_per_axis::Integer,
	output_start::Integer,
	digits::Integer,
	variable::AbstractString,
	fix_spec::AbstractString = "",
)
	incar = parse_incar(input_path)

	if :MAGMOM in keys(incar)
		magmom = incar[:MAGMOM]
	elseif :M_CONSTR in keys(incar)
		magmom = incar[:M_CONSTR]
	else
		error("MAGMOM or M_CONSTR not found in INCAR file")
	end

	if magmom isa AbstractString
		magmom_cleaned = replace(magmom, r"\\\s*[\r\n]+" => " ")
		magmom_cleaned = replace(magmom_cleaned, "\\" => "")
		magmom_cleaned = replace(magmom_cleaned, r"[\r\n]+" => " ")
		magmom_cleaned = strip(magmom_cleaned)
		magmom = [parse(Float64, n) for n in filter(!isempty, split(magmom_cleaned, r"\s+"))]
	elseif magmom isa AbstractVector
		magmom = [Float64(x) for x in magmom]
	else
		error("Invalid MAGMOM/M_CONSTR type: $(typeof(magmom))")
	end

	if length(magmom) % 3 != 0
		error("Invalid MAGMOM/M_CONSTR length: $(length(magmom)); must be a multiple of 3")
	end

	num_atoms = length(magmom) ÷ 3
	spin_matrix = zeros(3, num_atoms)
	for i in 1:num_atoms
		spin_matrix[:, i] = magmom[(3i-2):3i]
	end

	fixed_atom_indices = parse_atom_index_spec(fix_spec; max_index = num_atoms)
	n_axes = length(axes)

	for (axis_idx, axis) in enumerate(axes)
		for rep in 1:n_per_axis
			if variable == "tau"
				output_spin_matrix = mfa_spin_sampler_with_τ(spin_matrix, value)
			elseif variable == "m"
				output_spin_matrix = mfa_spin_sampler_with_magnetization(spin_matrix, value)
			else
				error("Invalid variable: $variable. Expected \"tau\" or \"m\".")
			end

			R = rotation_matrix_from_vectors([0.0, 0.0, 1.0], axis)
			output_spin_matrix = R * output_spin_matrix

			for idx in fixed_atom_indices
				@inbounds output_spin_matrix[:, idx] = R * spin_matrix[:, idx]
			end

			output_spin_flattened = reshape(output_spin_matrix, :)
			output_incar = copy(incar)
			output_incar[:MAGMOM] = output_spin_flattened
			output_incar[:M_CONSTR] = output_spin_flattened

			sample_id = output_start + (axis_idx - 1) * n_per_axis + rep
			output_path = @sprintf("sample-%0*d.INCAR", digits, sample_id)
			try
				write_incar(output_path, output_incar; wrap_vectors = true)
			catch e
				@error "Failed to write INCAR file" path = output_path exception = (e, catch_backtrace())
				rethrow(e)
			end
		end
	end
end

function main()
	s = ArgParseSettings(description = """
	MFA spin sampling with quantization axes along Lebedev grid directions.
	Specify a single value for tau or m. Grid: generate with --order or load from --lebedev file.
	""")

	@add_arg_table s begin
		"input"
		help = "Input INCAR path"
		required = true

		"variable"
		help = "Variable: tau or m"
		range_tester = in(["tau", "m"])
		required = true

		"value"
		help = "Single value (T/Tc for tau, or magnetization for m)"
		required = true
		arg_type = Float64

		"--order"
		help = "Lebedev order (generate grid in script; e.g. 9 -> 38 points). If set, --lebedev is ignored."
		arg_type = Int
		default = 0

		"--num_samples", "-n"
		help = "Number of samples per Lebedev direction (default: 1)"
		arg_type = Int
		default = 1

		"--lebedev", "-l"
		help = "Path to Lebedev grid file (used only when --order is 0; default: grid.txt)"
		arg_type = String
		default = "grid.txt"

		"--theta-min"
		help = "Polar angle lower bound [deg] (0 = +z); default: 0"
		arg_type = Float64
		default = 0.0

		"--theta-max"
		help = "Polar angle upper bound [deg]; default: 180"
		arg_type = Float64
		default = 180.0

		"--phi-min"
		help = "Azimuth lower bound [deg]; default: 0"
		arg_type = Float64
		default = 0.0

		"--phi-max"
		help = "Azimuth upper bound [deg]; default: 360"
		arg_type = Float64
		default = 360.0

		"--output-grid", "-o"
		help = "Write (filtered) grid points to FILE (x y z [weight]). No output if not set."
		arg_type = String
		default = ""

		"--fix"
		help = "Fix magnetic moments for 1-based atom indices (e.g. \"1-10,12\"). Same uniform rotation as other atoms."
		arg_type = String
		default = ""
	end

	args = parse_args(s)
	value = args["value"]
	n_per_axis = args["num_samples"]
	n_per_axis > 0 || error("--num_samples (-n) must be positive")

	if args["order"] > 0
		axes, weights = generate_lebedev_grid(args["order"])
		grid_src = "order $(args["order"]) ($(length(axes)) points)"
	else
		axes, weights = load_lebedev_axes(args["lebedev"])
		grid_src = args["lebedev"]
	end

	axes, weights = filter_axes_by_angle(
		axes, weights,
		Float64(args["theta-min"]), Float64(args["theta-max"]),
		Float64(args["phi-min"]), Float64(args["phi-max"]),
	)
	isempty(axes) && error("No grid points in the specified angular range")

	if !isempty(args["output-grid"])
		write_grid_file(args["output-grid"], axes, weights)
		@printf("Wrote %d grid points to %s\n", length(axes), args["output-grid"])
	end

	n_axes = length(axes)
	total_files = n_axes * n_per_axis
	digits = length(string(total_files))

	@printf("Lebedev grid: %s (%d directions after angle filter)\n", grid_src, n_axes)
	@printf("  theta: [%.1f, %.1f] deg  phi: [%.1f, %.1f] deg\n",
		args["theta-min"], args["theta-max"], args["phi-min"], args["phi-max"])
	@printf("Input: %s  Variable: %s = %.4f\n", args["input"], args["variable"], value)
	@printf("Samples per direction (-n): %d  Total INCAR files: %d\n", n_per_axis, total_files)

	sampling_mfa_lebedev(
		args["input"],
		value,
		axes,
		n_per_axis,
		0,
		digits,
		args["variable"],
		args["fix"],
	)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
