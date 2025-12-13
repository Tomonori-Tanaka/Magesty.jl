module EnergyTorque
using EzXML
using LinearAlgebra
using StaticArrays
using ..Symmetries
using ..Basis
using ..Optimize

export calc_energy, calc_torque

"""
	calc_energy(salc_list, spin_directions, symmetry, optimize)

Calculate the energy of a spin configuration using the spin cluster expansion.

# Description
Computes the energy for a given spin configuration by evaluating the spin cluster expansion
model. The energy is calculated as the dot product of the design vector (computed from the
spin directions and basis functions) with the optimized SCE coefficients, plus the reference energy.

# Arguments
- `salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}}`: List of symmetry-adapted linear combinations (SALCs), where each element is a vector of coupled basis functions belonging to the same key group.
- `spin_directions::AbstractMatrix{<:Real}`: A 3×N matrix representing the spin configuration, where each column is a unit vector (x, y, z) representing the spin direction at each atom. N is the number of atoms in the supercell.
- `symmetry::Symmetry`: Symmetry information of the crystal structure.
- `optimize::Optimizer`: Optimizer instance containing the optimized SCE coefficients (`optimize.SCE`) and reference energy (`optimize.reference_energy`).

# Returns
- `Float64`: The calculated energy of the spin configuration in eV.

# Throws
- `ArgumentError`: If `spin_directions` is not a 3×N matrix.

# Example
```julia
# Calculate energy for a spin configuration
spin_config = rand(3, num_atoms)
# Normalize to unit vectors
for i in 1:num_atoms
    spin_config[:, i] ./= norm(spin_config[:, i])
end
energy = calc_energy(salc_list, spin_config, symmetry, optimizer)
```
"""
function calc_energy(
    salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
    spin_directions::AbstractMatrix{<:Real},
    symmetry::Symmetry,
    optimize::Optimizer,
)::Float64
    if size(spin_directions, 1) != 3
        throw(ArgumentError("spin_directions must be a 3xN matrix"))
    end

	num_salcs = length(salc_list)
	design_vector = Vector{Float64}(undef, num_salcs)

	for i in 1:num_salcs
		key_group::Vector{Basis.CoupledBasis_with_coefficient} = salc_list[i]
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factor = (4*pi)^(n_C/2)  # (√(4π))^{n_C}
		
		# Sum contributions from all CoupledBasis_with_coefficient in this key group
		group_value = 0.0
		for cbc::Basis.CoupledBasis_with_coefficient in key_group
			group_value += Optimize.design_matrix_energy_element(
				cbc,
				spin_directions,
				symmetry,
			)
		end
		design_vector[i] = group_value * scaling_factor
	end

	return dot(design_vector, optimize.SCE) + optimize.reference_energy
end

"""
	calc_torque(salc_list, spin_directions, symmetry, optimize)

Calculate the torque for each atom in a spin configuration using the spin cluster expansion.

# Description
Computes the torque for each atom by evaluating the gradient of the spin cluster
expansion model with respect to the spin directions. The torque is calculated as the cross product of
the spin direction with the gradient of the energy with respect to the spin direction.

# Arguments
- `salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}}`: List of symmetry-adapted linear combinations (SALCs), where each element is a vector of coupled basis functions belonging to the same key group.
- `spin_directions::AbstractMatrix{<:Real}`: A 3×N matrix representing the spin configuration, where each column is a unit vector (x, y, z) representing the spin direction at each atom. N is the number of atoms in the supercell.
- `symmetry::Symmetry`: Symmetry information of the crystal structure.
- `optimize::Optimizer`: Optimizer instance containing the optimized SCE coefficients (`optimize.SCE`).

# Returns
- `Matrix{Float64}`: A 3×N matrix representing the torque for each atom, where each column corresponds to the torque vector (in eV) for one atom. The torque is the component of the local magnetic field perpendicular to the spin direction.

# Throws
- `ArgumentError`: If `spin_directions` is not a 3×N matrix.

# Example
```julia
# Calculate torque for a spin configuration
spin_config = rand(3, num_atoms)
# Normalize to unit vectors
for i in 1:num_atoms
    spin_config[:, i] ./= norm(spin_config[:, i])
end
torque = calc_torque(salc_list, spin_config, symmetry, optimizer)
# torque[:, i] gives the torque vector for atom i
```
"""
function calc_torque(
    salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
    spin_directions::AbstractMatrix{<:Real},
    symmetry::Symmetry,
    optimize::Optimizer,
)::Matrix{Float64}
    if size(spin_directions, 1) != 3
        throw(ArgumentError("spin_directions must be a 3xN matrix"))
    end

	num_atoms = size(spin_directions, 2)
	num_salcs = length(salc_list)
	torque = zeros(Float64, 3, num_atoms)

	@inbounds for iatom in 1:num_atoms
		@views dir_iatom = spin_directions[:, iatom]
		@inbounds for (salc_idx, key_group) in enumerate(salc_list)
			# Sum contributions from all CoupledBasis_with_coefficient in this key group
			group_grad = MVector{3, Float64}(0.0, 0.0, 0.0)
			for cbc in key_group
				grad_u = Optimize.calc_∇ₑu(cbc, iatom, spin_directions, symmetry)
				group_grad .+= grad_u
			end
			n_C = length(key_group[1].atoms)  # Number of sites in the cluster
			scaling_factor = (4*pi)^(n_C/2)  # (√(4π))^{n_C}
			torque[:, iatom] .+= cross(dir_iatom, Vector{Float64}(group_grad)) .* scaling_factor .* optimize.SCE[salc_idx]
		end
	end

	return torque
end


end
