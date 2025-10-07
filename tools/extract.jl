using ArgParse
using Printf
using Random

function extract_energy_magmom(
	file::AbstractString,
	energy_kind::AbstractString,
	mint::Bool,
)::Tuple{Float64, Matrix{Float64}}
	magmom_listoflist = Vector{Vector{Float64}}()
	magmom_temp_listoflist = Vector{Vector{Float64}}()

	energy::Float64 = 0.0
	energy_found = false
	open(file, "r") do f
		collecting = false
		for line in eachline(f)
			# extract magmom
			if occursin("M_int", line)
				collecting = true
				empty!(magmom_temp_listoflist)
				continue
			end

			if collecting && (occursin(":", line)) || (length(split(line)) ≠ 7) # end of magmom section (":" is intended to such as "DAV:" or "RMM:" keyword)
				collecting = false
				magmom_listoflist = deepcopy(magmom_temp_listoflist)
			end

			if collecting
				parts = split(line)
				if mint
					moment_x = parse(Float64, parts[5])
					moment_y = parse(Float64, parts[6])
					moment_z = parse(Float64, parts[7])
					push!(magmom_temp_listoflist, [moment_x, moment_y, moment_z])
				else
					moment_x = parse(Float64, parts[2])
					moment_y = parse(Float64, parts[3])
					moment_z = parse(Float64, parts[4])
					push!(magmom_temp_listoflist, [moment_x, moment_y, moment_z])
				end
			end

			# extract energy
			if occursin("F=", line)
				if energy_kind == "f"
					energy = parse(Float64, split(line)[3])
					energy_found = true
				elseif energy_kind == "e0"
					energy = parse(Float64, split(line)[5])
					energy_found = true
				end
			end
		end
	end
	if !energy_found
		error("energy not found in $file")
	end
	magmom_matrix::Matrix{Float64} = reduce(vcat, magmom_listoflist')
	return energy, magmom_matrix
end

function extract_magfield(file::String, num_atoms::Int)::Matrix{Float64}
	magfield_matrix = zeros(Float64, num_atoms, 3)

	open(file, "r") do f
		collecting = false
		atom_count = 0

		collecting = false
		for line in eachline(f)
			if occursin("lambda*MW_perp", line)
				collecting = true
				continue
			end

			if collecting
				parts = split(line)
				atom_idx = parse(Int, parts[1])
				magfield_x = parse(Float64, parts[2])
				magfield_y = parse(Float64, parts[3])
				magfield_z = parse(Float64, parts[4])
				magfield_matrix[atom_idx, :] = [magfield_x, magfield_y, magfield_z]
				atom_count += 1
			end

			if atom_count == num_atoms
				break
			end
		end
	end

	return magfield_matrix
end

function print_embset(
	concated_data::Matrix{Float64},
	energy::Float64,
	idx::Int,
	file::AbstractString,
)
	println(
		"# $idx, $file, energy unit = eV, magmom unit = Bohr magneton, magnetic field unit = T",
	)
	@printf("%.5f\n", energy)
	for (row_idx, row) in enumerate(eachrow(concated_data))
		println(
			lpad(row_idx, 7),
			"  ",
			join([lpad(@sprintf("%.7f", x), 12) for x in row[1:3]], " "),
			"  ",
			join([lpad(@sprintf("%.5e", x), 14) for x in row[4:6]], " "),
		)
	end
end


s = ArgParseSettings(
	description = """
Extract the energy, magnetic moments, and local magnetic field from the OUTCAR file and print the result in the EMBSET format.
""",
)
@add_arg_table s begin
	"--energy_kind", "-e"
	help = "The kind of energy. f: free energy, e0: energy(sigma->0)"
	default = "f"

	"--mint"
	help = "flag to extract magnetic moment from \"M_int\". If not specified, the magnetic moment is extracted from \"MW_int\". Usually, you should not use this flag."
	action = :store_true

	"--saxis", "-s"
	help = "quantization axis (three real numbers)"
	nargs = 3
	arg_type = Float64
	default = [0.0, 0.0, 1.0]

	"--target_files", "-f"
	help = "target files"
	required = true
	nargs = '+'

	"--randomize"
	help = "randomize the order of the spin configurations"
	action = :store_true
end

parsed_args = parse_args(ARGS, s)

# Randomize target files if --randomize flag is specified
target_files = parsed_args["target_files"]
if parsed_args["randomize"]
	shuffle!(target_files)
end

for (idx, file) in enumerate(target_files)
	energy, magmom_matrix =
		extract_energy_magmom(file, parsed_args["energy_kind"], parsed_args["mint"])
	num_atoms = size(magmom_matrix, 1)
	magfield_matrix = extract_magfield(file, num_atoms)

	# calculate rotation matrix
	alpha = atan(parsed_args["saxis"][2], parsed_args["saxis"][1])
	beta = atan(
		sqrt(parsed_args["saxis"][1]^2 + parsed_args["saxis"][2]^2),
		parsed_args["saxis"][3],
	)
	Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
	Ry(b) = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
	R = Rz(alpha) * Ry(beta)

	# rotate magmom_matrix and magfield_matrix
	magmom_matrix = (R * magmom_matrix')'
	magfield_matrix = (R * magfield_matrix')'

	concated_data = hcat(magmom_matrix, magfield_matrix)
	print_embset(concated_data, energy, idx, file)
end
