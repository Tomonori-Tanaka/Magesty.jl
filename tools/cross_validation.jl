module CrossValidation

using Random

include("../src/Magesty.jl")


"""

# Arguments
- `weight_min::Real`: Minimum weight 
- `weight_max::Real`: Maximum weight
- `weight_step::Integer`: Step size for weight
- `k_fold::Integer`: Number of folds for cross validation
- `initial_sce::AbstractVector{<:Real}`: Initial SCE coefficients
- `intitlal_ref_energy::Real`: Initial reference energy
- `shuffle::Bool`: Whether to shuffle the data set

# Returns
- `result_SCE::Vector{Float64}`: Best SCE coefficients to minimize the RMSE for test data
- `result_ref_energy::Float64`: Best reference energy to minimize the RMSE for test data
- `result_weight::Float64`: Best weight to minimize the RMSE for test data
- `result_RMSE4test_list::Vector{Float64}`: List of cv scores for test data
- `result_RMSE4train_list::Vector{Float64}`: List of cv scores for train data
"""

function cross_validation(
	optimize::SCEOptimizer,
	weight_min::Real,
	weight_max::Real,
	weight_step::Integer,
	k_fold::Integer,
	initial_sce::AbstractVector{<:Real},
	intitlal_ref_energy::Real,
	;
	shuffle_data::Bool = true,
)
	spinconfig_list::Vector{SpinConfig} = optimize.spinconfig_dataset.spinconfigs
	splitted_dataset::Vector{Vector{SpinConfig}} = split_data(spinconfig_list, k_fold; shuffle_data = shuffle_data)
	println("N: ", length(spinconfig_list), " k_fold: ", k_fold)
	for i in 1:k_fold
		println("Fold ", i, ": ", length(splitted_dataset[i]))
	end
end

function cross_validation_with_a_weight(sc::SpinCluster, weight::Real, splitted_dataset::Vector{Vector{SpinConfig}})
	for (i, test_list::Vector{SpinConfig}) in enumerate(splitted_dataset)
		train_list::Vector{SpinConfig} = vcat(splitted_dataset[1:i-1]..., splitted_dataset[i+1:end]...)
		energy_list_train = Vector{Float64}(undef, length(train_list))
		for (j, train_spinconfig::SpinConfig) in enumerate(train_list)
			sc_train 
			energy_list_train[j] = Magesty.calc_energy(sc, train_spinconfig.spin_directions)
		end
	end
end

"""
	split_data(spinconfig_list, k_fold; shuffle=true)

Split the dataset for k-fold cross validation.

# Arguments
- `spinconfig_list::AbstractVector{T}`: Dataset to be split
- `k_fold::Integer`: Number of folds
- `shuffle::Bool`: Whether to shuffle the data (default: true)

# Returns
- `Vector{Vector{T}}`: Dataset split into k subsets

# Throws
- `ArgumentError`: If k_fold is less than 2 or larger than the dataset size
"""
function split_data(
	spinconfig_list::AbstractVector{T},
	k_fold::Integer;
	shuffle_data::Bool = true,
) where T
	# Input validation
	if k_fold < 2
		throw(ArgumentError("k_fold must be at least 2"))
	end

	N = length(spinconfig_list)
	if k_fold > N
		throw(ArgumentError("k_fold ($k_fold) cannot be larger than the dataset size ($N)"))
	end

	# Split the data
	vec_indices::Vector{Int} = collect(1:N)
	if shuffle_data
		shuffle!(vec_indices)
	end
	base, rem = divrem(N, k_fold)
	indices = Vector{Vector{Int}}()
	idx = 1
	for i in 1:k_fold
		n = base + (i <= rem ? 1 : 0)
		push!(indices, vec_indices[idx:(idx+n-1)])
		idx += n
	end

	# Create split datasets
	return [spinconfig_list[list] for list in indices]::Vector{Vector{T}}
end

end
