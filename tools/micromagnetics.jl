"""
Calculation tool for the derivation of the micromagnetics model parameters.
"""

using ArgParse
using LinearAlgebra
using EzXML

using Magesty

@isdefined(convert2tensor) || begin
	include("./convert2tensor.jl")
	using .ExchangeTensor
end


function calc_micromagnetics(
	input_xml::String,
	basis::Magesty.SCEBasis,
	cutoff::Union{Float64, Nothing} = nothing,
)::Tuple{Matrix{Float64}, Matrix{Float64}}
	num_atoms = basis.structure.supercell.num_atoms
	atoms_in_prim = basis.symmetry.atoms_in_prim   # atom indices in the primitive cell
	min_distance_pairs = Magesty.Clusters._set_mindist_pairs(
		num_atoms,
		basis.structure.x_image_cart,
		basis.structure.exist_image,
	)

	stiffness_matrix = zeros(Float64, 3, 3)
	spiralization_matrix = zeros(Float64, 3, 3)

	for idx in 1:length(atoms_in_prim)
		i_atom = atoms_in_prim[idx]

		for i_pair in 1:num_atoms
			if i_pair == i_atom
				continue
			end

			# Get distance information (used for both cutoff check and calculation)
			dist_info_list = min_distance_pairs[i_atom, i_pair]

			if cutoff !== nothing
				distance = dist_info_list[1].distance
				if distance > cutoff
					continue
				end
			end

			# calculate the exchange interaction tensor
			exchange_tensor = convert2tensor(input_xml, [i_atom, i_pair])
			# Convert from <i, j> (per bond) to (i, j) (allow for double counting)
			exchange_tensor = (1 / 2) * exchange_tensor
			jij = exchange_tensor.isotropic_jij
			dm_vector = exchange_tensor.dm_vector
			for dist_info in dist_info_list
				# relative vector from i_atom to i_pair in Cartesian coordinates
				relvec::Vector{Float64} = dist_info.relative_vector

				# Avoid creating temporary arrays by using direct element access
				for i in 1:3, j in 1:3
					stiffness_matrix[i, j] += 0.5 * jij * relvec[i] * relvec[j]
				end

				for i in 1:3, j in 1:3
					spiralization_matrix[i, j] += dm_vector[i] * relvec[j]
				end
			end
		end
	end

	return stiffness_matrix, spiralization_matrix
end

function main()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--input_xml", "-x"
		help = "Input XML file with SCE basis and coefficients (jphi.xml)"
		required = true

		"--cutoff", "-c"
		help = "Cutoff distance for atom pairs (in Angstrom). If not specified, all atom pairs are considered"
		default = nothing
	end

	args = parse_args(s)
	basis = Magesty.load(Magesty.SCEBasis, args["input_xml"])

	cutoff = args["cutoff"] === nothing ? nothing : parse(Float64, args["cutoff"])
	stiffness_matrix, spiralization_matrix =
		calc_micromagnetics(args["input_xml"], basis, cutoff)
	println("Stiffness matrix:")
	display(stiffness_matrix)
	println("")
	println("Spiralization matrix:")
	display(spiralization_matrix)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
