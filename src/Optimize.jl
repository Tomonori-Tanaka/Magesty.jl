"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

using Base.Threads
using LinearAlgebra
using Printf
using MultivariateStats
using Statistics
using StaticArrays
using ..MySphericalHarmonics
using ..SALCs
using ..AtomicIndices
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..Basis
using ..BasisSets
using ..SpinConfigs

export Optimizer

struct Optimizer
	spinconfig_list::Vector{SpinConfig}
	reference_energy::Float64
	SCE::Vector{Float64}
	metrics::Dict{Symbol, Any}
	predicted_energy_list::Vector{Float64}
	predicted_torque_list::Vector{Matrix{Float64}}

	function Optimizer(
		structure::Structure,
		symmetry::Symmetry,
		basisset::BasisSet,
		alpha::Real,
		lambda::Real,
		weight::Real,
		spinconfig_list::AbstractVector{SpinConfig},
		;
		verbosity::Bool = true,
	)
		# Start timing
		start_time = time_ns()

		if verbosity
			println("""

			OPTIMIZATION
			============

			""")
		end

		# construct design matrix for energy and torque
		# salc_list is Vector{Vector{...}} where each inner vector represents one key group
		# design_matrix should have width = length(salc_list) + 1 (one column per key group + bias)
		if verbosity
			println("Constructing design matrix for energy...")
		end
		design_matrix_energy = build_design_matrix_energy(
			basisset.salc_list,
			spinconfig_list,
			symmetry,
		)
		if verbosity
			println("Constructing design matrix for torque...")
		end

		design_matrix_torque = build_design_matrix_torque(
			basisset.salc_list,
			spinconfig_list,
			structure.supercell.num_atoms,
			symmetry,
		)

		# construct observed_energy_list and observed_torque_list
		if verbosity
			println("Constructing observed data for energy and torque...")
		end
		observed_energy_list = [spinconfig.energy for spinconfig in spinconfig_list]
		observed_torque_list = [spinconfig.torques for spinconfig in spinconfig_list]

		if verbosity
			println("Fitting SCE coefficients...\n")
		end
		j0, jphi = elastic_net_regression(
			design_matrix_energy,
			design_matrix_torque,
			observed_energy_list,
			observed_torque_list,
			alpha,
			lambda,
			weight,
		)

		predicted_energy_list = design_matrix_energy[:, 2:end] * jphi .+ j0
		predicted_torque_flattened_list::Vector{Float64} = design_matrix_torque * jphi

		# Reshape predicted_torque_flattened_list to a vector of matrices
		block_size = 3 * structure.supercell.num_atoms
		num_configs = length(spinconfig_list)
		predicted_torque_list::Vector{Matrix{Float64}} = [
			reshape(
				predicted_torque_flattened_list[((i-1)*block_size+1):(i*block_size)],
				3,
				structure.supercell.num_atoms,
			)
			for i in 1:num_configs
		]

		metrics = calc_metrics(
			observed_energy_list,
			predicted_energy_list,
			observed_torque_list,
			predicted_torque_list,
		)

		if verbosity
			print_sce_coeffs(j0, jphi)
			print_metrics(metrics)

			println(@sprintf(
				"""

				Time Elapsed: %.6f sec.
				""",
				(time_ns() - start_time) / 1e9,
			))
		end

		return new(
			spinconfig_list,
			j0,
			jphi,
			metrics,
			predicted_energy_list,
			predicted_torque_list,
		)
	end
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	alpha::Real,
	lambda::Real,
	weight::Real,
	datafile::AbstractString,
	;
	verbosity::Bool = true,
)
	# read datafile
	spinconfig_list = read_embset(datafile, structure.supercell.num_atoms)

	return Optimizer(
		structure,
		symmetry,
		basisset,
		alpha,
		lambda,
		weight,
		spinconfig_list,
		verbosity = verbosity,
	)
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	config::Config4Optimize,
	;
	verbosity::Bool = true,
)
	return Optimizer(
		structure,
		symmetry,
		basisset,
		config.alpha,
		config.lambda,
		config.weight,
		config.datafile,
		verbosity = verbosity,
	)
end


function build_design_matrix_energy(
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs + 1)

	# set first column to 1 (reference_energy term)
	design_matrix[:, 1] .= 1.0

	for i in 1:num_salcs
		key_group::Vector{Basis.CoupledBasis_with_coefficient} = salc_list[i]
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factor = (4*pi)^(n_C/2)  # (√(4π))^{n_C}
		@inbounds for j in 1:num_spinconfigs
			# Sum contributions from all CoupledBasis_with_coefficient in this key group
			group_value = 0.0
			for cbc::Basis.CoupledBasis_with_coefficient in key_group
				group_value += design_matrix_energy_element(
					cbc,
					spinconfig_list[j].spin_directions,
					symmetry,
				)
			end
			design_matrix[j, i+1] = group_value * scaling_factor
		end
	end

	return design_matrix
end


"""
	build_design_matrix_energy(salc_list, spinconfig_list, symmetry) -> Matrix{Float64}

Build the energy design matrix used for regression.

# Description
- Returns a matrix of size (num_configs × (num_salcs + 1)).
- The first column is the bias term (all ones). Columns 2..end are features per SALC.

# Arguments
- `salc_list::AbstractVector{SALC}`: List of SALC objects
- `spinconfig_list::AbstractVector{SpinConfig}`: Vector of SpinConfig objects
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Matrix{Float64}`: Energy design matrix
"""
function build_design_matrix_energy(
	salc_list::AbstractVector{SALC},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs + 1)

	# set first column to 1 (reference_energy term)
	design_matrix[:, 1] .= 1.0

	for i in 1:num_salcs
		@inbounds for j in 1:num_spinconfigs
			design_matrix[j, i+1] = design_matrix_energy_element(
				salc_list[i],
				spinconfig_list[j].spin_directions,
				symmetry,
			)
		end
	end

	return design_matrix
end

"""
	design_matrix_energy_element(salc, spin_directions, symmetry) -> Float64

Compute one energy-design feature for a given SALC and spin directions.

# Description
- Contracts spherical harmonics over atoms following symmetry translations.
- Equivalent to one column entry (excluding bias) in the energy design matrix.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure

# Returns
- `Float64`: Feature value for the SALC
"""
function design_matrix_energy_element(
	salc::SALC,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Float64

	result::Float64 = 0.0
	for itrans in symmetry.symnum_translation
		@inbounds for (basis_idx, basis::IndicesUniqueList) in enumerate(salc.basisset)
			product_tmp::Float64 = 1.0
			@inbounds for ibasis::Indices in basis
				atom::Int = symmetry.map_sym[ibasis.atom, itrans]
				l::Int = ibasis.l
				m::Int = ibasis.m
				product_tmp *= @views Sₗₘ(l, m, spin_directions[:, atom])
			end
			result += salc.coeffs[basis_idx] * salc.multiplicity[basis_idx] * product_tmp
		end
	end

	return result
end

"""
	design_matrix_energy_element(cbc, spin_directions, symmetry) -> Float64

Compute one energy-design feature for a given CoupledBasis_with_coefficient and spin directions.

# Description
- Contracts coupled angular momentum basis tensor with spherical harmonics over atoms following symmetry translations.
- Equivalent to one column entry (excluding bias) in the energy design matrix.

# Arguments
- `cbc::Basis.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure

# Returns
- `Float64`: Feature value for the CoupledBasis_with_coefficient
"""
function design_matrix_energy_element(
	cbc::Basis.CoupledBasis_with_coefficient,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Float64
	result::Float64 = 0.0
	N = length(cbc.atoms)

	for itrans in symmetry.symnum_translation
		# Translate atoms
		translated_atoms = [symmetry.map_sym[atom, itrans] for atom in cbc.atoms]

		# Compute spherical harmonics for each site
		sh_values = Vector{Vector{Float64}}(undef, N)
		for (site_idx, atom) in enumerate(translated_atoms)
			l = cbc.ls[site_idx]
			sh_values[site_idx] = Vector{Float64}(undef, 2*l+1)
			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1
				sh_values[site_idx][m_idx] = @views Sₗₘ(l, m, spin_directions[:, atom])
			end
		end

		# Contract coeff_tensor with spherical harmonics
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size)
		# where di = 2*li + 1
		tensor_result = 0.0
		Mf_size = size(cbc.coeff_tensor, N+1)
		dims = [2*l + 1 for l in cbc.ls]

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_contribution = 0.0

			# Iterate over all combinations of m indices using CartesianIndices
			# Create indices for first N dimensions
			site_indices = CartesianIndices(Tuple(dims))
			for site_idx_tuple in site_indices
				# Compute product of spherical harmonics
				product = 1.0
				for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
					product *= sh_values[site_idx][m_idx]
				end

				# Access tensor element: coeff_tensor[site_idx_tuple..., mf_idx]
				tensor_idx = (site_idx_tuple.I..., mf_idx)
				mf_contribution += cbc.coeff_tensor[tensor_idx...] * product
			end

			tensor_result += cbc.coefficient[mf_idx] * mf_contribution
		end

		result += tensor_result * cbc.multiplicity
	end

	return result
end

"""
	build_design_matrix_torque(salc_list, spinconfig_list, num_atoms, symmetry) -> Matrix{Float64}

Build the torque design matrix used for regression.

# Description
- For each spin configuration, constructs a block of size (3·num_atoms × num_salcs)
  whose rows are per-atom XYZ components of `cross(spin_dir, ∇ₑu)`.
- Blocks are vertically concatenated across configurations.

# Arguments
- `salc_list`: List of SALCs or CoupledBasis_with_coefficient
- `spinconfig_list`: Vector of spin configurations
- `num_atoms`: Number of atoms in the structure
- `symmetry`: Symmetry information

# Returns
- `Matrix{Float64}`: Torque design matrix
"""
function build_design_matrix_torque(
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	num_atoms::Integer,
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)

	design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)

		for (sc_idx, spinconfig) in enumerate(spinconfig_list)
		torque_design_block = zeros(Float64, 3*num_atoms, num_salcs)
		@inbounds for iatom in 1:num_atoms
			@views dir_iatom = spinconfig.spin_directions[:, iatom]
			@inbounds for (salc_idx, key_group) in enumerate(salc_list)
				# Sum contributions from all CoupledBasis_with_coefficient in this key group
				group_grad = MVector{3, Float64}(0.0, 0.0, 0.0)
				for cbc in key_group
					grad_u = calc_∇ₑu(cbc, iatom, spinconfig.spin_directions, symmetry)
					group_grad .+= grad_u
				end
				n_C = length(key_group[1].atoms)  # Number of sites in the cluster
				scaling_factor = (4*pi)^(n_C/2)  # (√(4π))^{n_C}
				@views torque_design_block[(3*(iatom-1)+1):(3*iatom), salc_idx] =
					cross(dir_iatom, Vector{Float64}(group_grad)) * scaling_factor
			end
		end
		design_matrix_list[sc_idx] = torque_design_block
	end

	return vcat(design_matrix_list...)
end

function build_design_matrix_torque(
	salc_list,
	spinconfig_list,
	num_atoms,
	symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)

	for (sc_idx, spinconfig) in enumerate(spinconfig_list)
		torque_design_block = zeros(Float64, 3*num_atoms, num_salcs)
		@inbounds for iatom in 1:num_atoms
			@views dir_iatom = spinconfig.spin_directions[:, iatom]
			@inbounds for (salc_idx, salc) in enumerate(salc_list)
				grad_u = calc_∇ₑu(salc, iatom, spinconfig.spin_directions, symmetry)
				@views torque_design_block[(3*(iatom-1)+1):(3*iatom), salc_idx] =
					cross(dir_iatom, grad_u)
			end
		end
		design_matrix_list[sc_idx] = torque_design_block
	end

	return vcat(design_matrix_list...)
end

"""
	calc_∇ₑu(salc, atom, spin_directions, symmetry) -> Vector{Float64}

Compute the gradient of the real spherical harmonic product for a SALC with respect to
the spin direction of a specific atom.

# Description
- Returns a 3-vector corresponding to (∂/∂x, ∂/∂y, ∂/∂z) components.
- Applies symmetry translations before accumulation.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination
- `atom::Integer`: Target atom index (1-based)
- `spin_directions::AbstractMatrix{<:Real}`: 3×N spin directions
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Vector{Float64}`: Gradient vector (length 3)
"""
function calc_∇ₑu(
	salc::SALC,
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Vector{Float64}
	spin_dirs = [SVector{3, Float64}(spin_directions[:, i]) for i in axes(spin_directions, 2)]

	result = MVector{3, Float64}(0.0, 0.0, 0.0)

	@inbounds for itrans in symmetry.symnum_translation
		translated_salc = translate_atom_idx_of_salc(salc, symmetry.map_sym, itrans)
		@inbounds for (basis, coeff, multiplicity) in
					  zip(
			translated_salc.basisset,
			translated_salc.coeffs,
			translated_salc.multiplicity,
		)
			# Skip if atom is not in the basis
			!any(idx -> idx.atom == atom, basis) && continue

			# Calculate product of spherical harmonics and their derivatives
			product = MVector{3, Float64}(1.0, 1.0, 1.0)
			@inbounds for indices in basis
				if indices.atom == atom
					product .*= ∂ᵢSlm(indices.l, indices.m, @inbounds spin_dirs[indices.atom])
				else
					product .*= Sₗₘ(indices.l, indices.m, @inbounds spin_dirs[indices.atom])
				end
			end
			result .+= product .* (coeff * multiplicity)
		end
	end
	return Vector{Float64}(result)
end

"""
	calc_∇ₑu(cbc, atom, spin_directions, symmetry) -> Vector{Float64}

Compute the gradient of the coupled angular momentum basis for a CoupledBasis_with_coefficient
with respect to the spin direction of a specific atom.

# Description
- Returns a 3-vector corresponding to (∂/∂x, ∂/∂y, ∂/∂z) components.
- Applies symmetry translations before accumulation.

# Arguments
- `cbc::Basis.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `atom::Integer`: Target atom index (1-based)
- `spin_directions::AbstractMatrix{<:Real}`: 3×N spin directions
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Vector{Float64}`: Gradient vector (length 3)
"""
function calc_∇ₑu(
	cbc::Basis.CoupledBasis_with_coefficient,
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Vector{Float64}
	result = MVector{3, Float64}(0.0, 0.0, 0.0)
	N = length(cbc.atoms)

	@inbounds for itrans in symmetry.symnum_translation
		# Translate atoms
		translated_atoms = [symmetry.map_sym[a, itrans] for a in cbc.atoms]

		# Check if atom is in the translated atoms list
		atom_site_idx = findfirst(==(atom), translated_atoms)
		if atom_site_idx === nothing
			continue
		end

		# Compute spherical harmonics and their derivatives for each site
		sh_values = Vector{Vector{Float64}}(undef, N)
		sh_grad_values = Vector{Vector{Vector{Float64}}}(undef, N)

		for (site_idx, translated_atom) in enumerate(translated_atoms)
			l = cbc.ls[site_idx]
			sh_values[site_idx] = Vector{Float64}(undef, 2*l+1)
			sh_grad_values[site_idx] = Vector{Vector{Float64}}(undef, 2*l+1)

			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1

				if site_idx == atom_site_idx
					# Use gradient for the target atom
					sh_grad_values[site_idx][m_idx] =
						@views ∂ᵢSlm(l, m, spin_directions[:, translated_atom])
					sh_values[site_idx][m_idx] =
						@views Sₗₘ(l, m, spin_directions[:, translated_atom])
				else
					# Use regular spherical harmonic for other atoms
					sh_values[site_idx][m_idx] =
						@views Sₗₘ(l, m, spin_directions[:, translated_atom])
					sh_grad_values[site_idx][m_idx] = [0.0, 0.0, 0.0]  # Not used, but needed for indexing
				end
			end
		end

		# Contract coeff_tensor with spherical harmonics and gradients
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size)
		# where di = 2*li + 1
		grad_result = MVector{3, Float64}(0.0, 0.0, 0.0)
		Mf_size = size(cbc.coeff_tensor, N+1)
		dims = [2*l + 1 for l in cbc.ls]

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_grad_contribution = MVector{3, Float64}(0.0, 0.0, 0.0)

			# Iterate over all combinations of m indices using CartesianIndices
			site_indices = CartesianIndices(Tuple(dims))
			for site_idx_tuple in site_indices
				# Compute product of spherical harmonics
				product = 1.0
				for (site_idx, m_idx) in enumerate(site_idx_tuple.I)
					if site_idx == atom_site_idx
						# Skip this site in the product (will multiply gradient separately)
						continue
					end
					product *= sh_values[site_idx][m_idx]
				end

				# Multiply by gradient for the target atom site
				m_idx_atom = site_idx_tuple.I[atom_site_idx]
				grad_atom = sh_grad_values[atom_site_idx][m_idx_atom]

				# Access tensor element: coeff_tensor[site_idx_tuple..., mf_idx]
				tensor_idx = (site_idx_tuple.I..., mf_idx)
				coeff_val = cbc.coeff_tensor[tensor_idx...]

				mf_grad_contribution .+= coeff_val * product .* grad_atom
			end

			grad_result .+= cbc.coefficient[mf_idx] .* mf_grad_contribution
		end

		result .+= grad_result .* cbc.multiplicity
	end

	return Vector{Float64}(result)
end

"""
	elastic_net_regression(design_matrix_energy, design_matrix_torque, observed_energy_list, observed_torque_list, alpha, lambda, weight)

Solve a combined regression for energy and torque using ridge (elastic net with alpha=0) regularization.

# Description
- Scales energy/torque parts by √weight and concatenates them into one system.
- Excludes the bias term from regularization.

# Arguments
- `design_matrix_energy`: Energy design matrix (bias column included)
- `design_matrix_torque`: Torque design matrix (no bias column)
- `observed_energy_list`: Observed energies
- `observed_torque_list`: Observed torques as 3×num_atoms matrices per configuration
- `alpha`: Unused (kept for API compatibility)
- `lambda`: Regularization strength
- `weight`: Trade-off between energy (1-weight) and torque (weight)

# Returns
- `(j0::Float64, jphi::Vector{Float64})`: Bias and coefficients
"""
function elastic_net_regression(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_torque::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
	alpha::Real,
	lambda::Real,
	weight::Real,
)
	# weight parameters
	w_e = 1 - weight
	w_m = weight

	# Flatten observed torque
	observed_torque_flattened = vcat(vec.(observed_torque_list)...)

	# Normalize the design matrices by using factor of 1/2N_data and √weight
	normalized_design_matrix_energy =
		design_matrix_energy .* sqrt(w_e)
	normalized_design_matrix_energy[:, 1] .= 1.0
	normalized_design_matrix_torque =
		design_matrix_torque .* sqrt(w_m)

	# Also normalise the observed vectors
	normalized_observed_energy_list =
		observed_energy_list .* sqrt(w_e)
	normalized_observed_torque_flattened =
		observed_torque_flattened .* sqrt(w_m)

	# Add 0 bias term to the design matrix for torque
	# to align with the energy design matrix
	with_bias_design_matrix_torque =
		hcat(zeros(size(normalized_design_matrix_torque, 1)), normalized_design_matrix_torque)

	# Construct the augmented design matrix
	X = vcat(
		normalized_design_matrix_energy,
		with_bias_design_matrix_torque,
	)

	y = vcat(
		normalized_observed_energy_list,
		normalized_observed_torque_flattened,
	)

	# betas = X \ y
	# jphi = betas[2:end]
	# j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	# println("j0: $j0")
	# println("jphi: $jphi")

	# Elastic net regression solution
	lambda_vec = fill(lambda, size(X, 2))
	lambda_vec[1] = 0.0  # exclude bias term from regularization
	j_values = begin
		if lambda ≈ 0.0
			X \ y
		else
			ridge(X, y, lambda_vec; bias = false)
		end
	end
	jphi = j_values[2:end]
	j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)


	# fit = glmnet(X, y; alpha = alpha, lambda = [lambda], standardize = true)
	# fit = glmnet(normalized_design_matrix_energy, normalized_observed_energy_list; alpha = alpha, lambda = [lambda], standardize = true)
	# Extract coefficients
	# jphi = fit.betas[2:end, 1]
	# j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	# j0 = fit.a0[1] * sqrt(2 * num_data) * num_atoms / sqrt(w_e)


	return j0, jphi
end


function translate_atom_idx_of_salc(
	salc::SALC,
	map_sym::AbstractMatrix{<:Integer},
	itrans::Integer,
)::SALC
	translated_basisset = Vector{IndicesUniqueList}()
	for basis::IndicesUniqueList in salc.basisset
		translated_basis = IndicesUniqueList()
		for indices::Indices in basis
			translated_atom_idx = map_sym[indices.atom, itrans]
			push!(
				translated_basis,
				Indices(translated_atom_idx, indices.l, indices.m, indices.cell),
			)
		end
		push!(translated_basisset, translated_basis)
	end

	return SALC(translated_basisset, salc.coeffs, salc.multiplicity)
end

function calc_rmse(list1::AbstractVector{<:Real}, list2::AbstractVector{<:Real})::Float64
	# Calculate the Root Mean Square Error (RMSE) between two lists
	if length(list1) != length(list2)
		throw(ArgumentError("The lengths of the two lists must be equal."))
	end
	return sqrt(mean((list1 .- list2) .^ 2))
end

function calc_r2score(
	observed_list::AbstractVector{<:Real},
	predicted_list::AbstractVector{<:Real},
)::Float64
	# R² = 1 - (SS_res / SS_tot)
	# SS_res = Σ(y_observed - y_predicted)²
	# SS_tot = Σ(y_observed - y_mean)²
	ss_res = sum((observed_list .- predicted_list) .^ 2)
	ss_tot = sum((observed_list .- mean(observed_list)) .^ 2)
	return 1 - ss_res / ss_tot
end

"""
	calc_metrics(observed_energy_list, predicted_energy_list, observed_torque_list, predicted_torque_list) -> Dict{Symbol,Any}

Compute RMSE and R² metrics for energy and torque.

# Returns
- `Dict{Symbol,Any}` with keys: `:rmse_energy`, `:rmse_torque`, `:r2score_energy`, `:r2score_torque`
"""
function calc_metrics(
	observed_energy_list::AbstractVector{<:Real},
	predicted_energy_list::AbstractVector{<:Real},
	observed_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
	predicted_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
)::Dict{Symbol, Any}
	observed_torque_flattened_list = vcat(vec.(observed_torque_list)...)
	predicted_torque_flattened_list = vcat(vec.(predicted_torque_list)...)
	return Dict(
		:rmse_energy => calc_rmse(observed_energy_list, predicted_energy_list),
		:rmse_torque => calc_rmse(observed_torque_flattened_list, predicted_torque_flattened_list),
		:r2score_energy => calc_r2score(observed_energy_list, predicted_energy_list),
		:r2score_torque =>
			calc_r2score(observed_torque_flattened_list, predicted_torque_flattened_list),
	)
end

function print_sce_coeffs(reference_energy::Float64, sce_coeffs::Vector{Float64})
	ndigit = ndigits(length(sce_coeffs))
	println("	SCE coefficients:")
	println(@sprintf("	  E_ref: % .10f eV", reference_energy))
	for (i, sce_coeff) in enumerate(sce_coeffs)
		println(@sprintf("    %*d: % .10f", ndigit, i, sce_coeff))
	end
	println()
end

function print_metrics(
	metrics::Dict{Symbol, Any},
)
	println(
		@sprintf(
			"""
			Root Mean Square Error (RMSE)
			-----------------------------
			RMSE for energy: %.6f meV
			RMSE for magnetic field: %.6f meV

			R^2 Score
			---------
			R^2 for energy: %.6f
			R^2 for magnetic field: %.6f
			""",
			metrics[:rmse_energy] * 1000,
			metrics[:rmse_torque] * 1000,
			metrics[:r2score_energy],
			metrics[:r2score_torque],
		)
	)
end

end

