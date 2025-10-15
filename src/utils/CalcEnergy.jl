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
		design_list[i] = Optimize.design_matrix_energy_element(
			salc_list[i],
			spin_config,
			symmetry,
		)
	end

	return dot(design_list, optimize.SCE) + optimize.reference_energy
end

function calc_torque(
	salc_list::AbstractVector{SALC},
	spin_config::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	optimize::Optimizer,
)::Matrix{Float64}
	if size(spin_config, 1) != 3
		throw(ArgumentError("spin_config must be a 3xN matrix"))
	end

	num_atoms = size(spin_config, 2)

	# 3Natom column vector
	torque_flattened = Optimize.build_design_matrix_torque(salc_list, spin_config, symmetry,) * optimize.SCE
	torque = reshape(torque_flattened, 3, num_atoms)

	return torque
end


end
