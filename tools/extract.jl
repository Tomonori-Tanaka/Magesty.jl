using ArgParse
using Printf

function extract_energy_magmom(
	file::String,
	energy_kind::String,
	mint::Bool,
)::Tuple{Float64, Matrix{Float64}}
	magmom_listoflist = Vector{Vector{Float64}}()
	magmom_temp_listoflist = Vector{Vector{Float64}}()
	energy::Float64 = 0.0

	open(file, "r") do f
		collecting = false
		for line in eachline(f)
			# extract magmom
			if occursin("M_int", line)
				collecting = true
				empty!(magmom_temp_listoflist)
				continue
			end

			if collecting && occursin("DAV:", line)
				collecting = false
				magmom_listoflist = deepcopy(magmom_temp_listoflist)
				continue
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
				elseif energy_kind == "e0"
					energy = parse(Float64, split(line)[5])
				end
			end
		end
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


s = ArgParseSettings()
@add_arg_table s begin
	"--energy_kind"
	help = "The kind of energy. f: free energy, e0: energy(sigma->0)"
	default = "f"

	"--mint"
	help = "flag to extract magnetic moment from \"M_int\""
	action = :store_true

	"--target_files", "-t"
	help = "target files"
	required = true
	nargs = '+'
end

parsed_args = parse_args(ARGS, s)

for (idx, file) in enumerate(parsed_args["target_files"])
	energy, magmom_matrix =
		extract_energy_magmom(file, parsed_args["energy_kind"], parsed_args["mint"])
	num_atoms = size(magmom_matrix, 1)
	magfield_matrix = extract_magfield(file, num_atoms)

	concated_data = hcat(magmom_matrix, magfield_matrix)
	println(
		"# $idx, $file, energy unit = eV, magmom unit = Bohr magneton, magnetic field unit = T",
	)
	println(energy)
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
