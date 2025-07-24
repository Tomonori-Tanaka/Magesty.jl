module CalcEnergy
using EzXML
using LinearAlgebra
using ..SALCs
using ..Symmetries
using ..Optimize

export calc_energy

function calc_energy(
	salc_list::AbstractVector{SALC},
	spin_config::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	optimize::Optimizer,
)::Float64
	if size(spin_config, 1) != 3
		throw(ArgumentError("spin_config must be a 3xN matrix"))
	end

	design_list = Vector{Float64}(undef, length(salc_list))

	for i in eachindex(salc_list)
		design_list[i] = Optimize.calc_X_element_energy(
			salc_list[i],
			spin_config,
			symmetry,
		)
	end

	return dot(design_list, optimize.SCE) + optimize.reference_energy
end

function calc_magfield(
	salc_list::AbstractVector{SALC},
	spin_config::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	optimize::Optimizer,
)::Matrix{Float64}
	if size(spin_config, 1) != 3
		throw(ArgumentError("spin_config must be a 3xN matrix"))
	end

	num_atoms = size(spin_config, 2)
	
	# Calculate magnetic field for each atom in each direction (x, y, z)
	magfield = zeros(3, num_atoms)
	
	for atom_idx in 1:num_atoms
		for direction in 1:3
			row_idx = 3 * (atom_idx - 1) + direction
			
			# Use calc_X_element_magfield for each SALC
			for (i, salc) in enumerate(salc_list)
				magfield[direction, atom_idx] += 
					Optimize.calc_X_element_magfield(salc, spin_config, symmetry, row_idx) * 
					optimize.SCE[i]
			end
		end
	end

	return magfield
end


end
