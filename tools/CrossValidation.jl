module CrossValidation

using Random

using Statistics
using ..Magesty.SpinConfigs
using ..Magesty

export cross_validation, loocv

"""
    loocv(spincluster, weights)

Perform Leave-One-Out Cross Validation (LOOCV) for the given spin cluster and weights.

# Arguments
- `spincluster::SpinCluster`: The spin cluster object containing the dataset and optimization settings.
- `weights::AbstractVector{<:Real}`: A vector of regularization weights to be evaluated.

# Returns
- `Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, SpinCluster}`: 
  A tuple containing the optimal weight, the list of weights, the test RMSE list, the train RMSE list, and the optimal SpinCluster.

# Throws
- `ArgumentError`: If `weights` is empty or contains invalid values.
"""
function loocv(
	spincluster::SpinCluster,
	weights::AbstractVector{<:Real},
)
	spinconfig_num::Int = length(spincluster.optimize.spinconfig_list)

	return cross_validation(spincluster, weights, spinconfig_num; shuffle_data = false)
end


"""
    cross_validation(spincluster, weights, k_fold; shuffle_data=true)

Perform k-fold cross validation for the given spin cluster and a set of regularization weights.

# Arguments
- `spincluster::SpinCluster`: The spin cluster object containing the dataset and optimization settings.
- `weights::AbstractVector{<:Real}`: A vector of regularization weights to be evaluated.
- `k_fold::Integer`: The number of folds for cross validation.
- `shuffle_data::Bool`: Whether to shuffle the data before splitting into folds (default: true).

# Returns
- `Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, SpinCluster}`: 
  A tuple containing the optimal weight, the list of weights, the test RMSE list, the train RMSE list, and the optimal SpinCluster.

# Throws
- `ArgumentError`: If `weights` is empty or contains invalid values.
- `ArgumentError`: If `k_fold` is less than 2 or greater than the dataset size.
"""
function cross_validation(
	spincluster::SpinCluster,
	weights::AbstractVector{<:Real},
	k_fold::Integer,
	;
	shuffle_data::Bool = true,
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, SpinCluster}
	weight_list::Vector{Float64} = collect(weights)
	if isempty(weight_list)
		throw(ArgumentError("weights must not be empty: $weights"))
	elseif !all(weight_list .>= 0.0)
		throw(ArgumentError("weights must be positive: $weights"))
	elseif !all(weight_list .<= 1.0)
		throw(ArgumentError("weights must be less than or equal to 1: $weights"))
	end

	if k_fold < 2
		throw(ArgumentError("k_fold ($k_fold) must be at least 2"))
	end

	spinconfig_list::Vector{SpinConfig} = spincluster.optimize.spinconfig_list
	if k_fold > length(spinconfig_list)
		throw(
			ArgumentError(
				"k_fold ($k_fold) cannot be larger than the dataset size ($(length(spinconfig_list)))",
			),
		)
	end

	splitted_spinconfig_list::Vector{Vector{SpinConfig}} =
		split_data(spinconfig_list, k_fold; shuffle_data = shuffle_data)

	current_spincluster = spincluster
	test_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(weight_list))
	train_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(weight_list))
	spincluster_list::Vector{SpinCluster} = Vector{SpinCluster}(undef, length(weight_list))
	for (i, weight) in enumerate(weight_list)
		cv_score_test::Float64, cv_score_train::Float64, current_spincluster::SpinCluster =
			cross_validation_with_a_weight(
				current_spincluster,
				weight,
				splitted_spinconfig_list,
			)
		test_rmse_list[i] = cv_score_test
		train_rmse_list[i] = cv_score_train
		spincluster_list[i] = current_spincluster
	end

	optimum_weight::Float64 = weight_list[argmin(test_rmse_list)]
	optimum_spincluster::SpinCluster = spincluster_list[argmin(test_rmse_list)]
	return optimum_weight, weight_list, test_rmse_list, train_rmse_list, optimum_spincluster
end

"""
    cross_validation_with_a_weight(spincluster, weight, splitted_dataset)

Perform cross validation for a single regularization weight over the provided data splits.

# Arguments
- `spincluster::SpinCluster`: The base spin cluster object.
- `weight::Real`: The regularization weight to be evaluated.
- `splitted_dataset::Vector{Vector{SpinConfig}}`: The dataset split into k folds (as produced by `split_data`).

# Returns
- `Tuple{Float64, Float64, SpinCluster}`: 
  A tuple containing the average test RMSE, the average train RMSE, and the SpinCluster trained on all data with the given weight.

# Throws
- No explicit exceptions are thrown by this function, but errors may propagate from underlying methods if input data is invalid.
"""
function cross_validation_with_a_weight(
	spincluster::SpinCluster,
	weight::Real,
	splitted_dataset::Vector{Vector{SpinConfig}},
)::Tuple{Float64, Float64, SpinCluster}
	test_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(splitted_dataset))
	train_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(splitted_dataset))

	for (i, test_spinconfigs::Vector{SpinConfig}) in enumerate(splitted_dataset)
		# Organize training data
		train_spinconfigs::Vector{SpinConfig} =
			vcat(splitted_dataset[1:(i-1)]..., splitted_dataset[(i+1):end]...)

		# Create spin cluster for training data
		trained_spin_cluster::SpinCluster =
			SpinCluster(spincluster, weight, train_spinconfigs, false)
		# Calculate RMSE for training data
		predicted_energies_train::Vector{Float64} =
			trained_spin_cluster.optimize.predicted_energy_list
		observed_energies_train::Vector{Float64} =
			[train_spinconfig.energy for train_spinconfig in train_spinconfigs]
		train_rmse_list[i] =
			sqrt(mean((predicted_energies_train .- observed_energies_train) .^ 2))

		# Calculate RMSE for test data
		predicted_energies_test::Vector{Float64} =
			Vector{Float64}(undef, length(test_spinconfigs))
		observed_energies_test::Vector{Float64} =
			[test_spinconfig.energy for test_spinconfig in test_spinconfigs]
		for (j, test_spinconfig::SpinConfig) in enumerate(test_spinconfigs)
			predicted_energies_test[j] =
				Magesty.calc_energy(trained_spin_cluster, test_spinconfig.spin_directions)
		end
		test_rmse_list[i] =
			sqrt(mean((predicted_energies_test .- observed_energies_test) .^ 2))
	end

	# Calculate CV scores for test and training data
	cv_score_test::Float64 = mean(test_rmse_list)
	cv_score_train::Float64 = mean(train_rmse_list)

	# Calculate the spin cluster with the weight using all data
	all_spinconfigs::Vector{SpinConfig} = vcat(splitted_dataset...)
	spincluster_with_a_weight::SpinCluster =
		SpinCluster(spincluster, weight, all_spinconfigs, false)
	return cv_score_test, cv_score_train, spincluster_with_a_weight
end

"""
	split_data(spinconfig_list, k_fold; shuffle_data=true)

Split the dataset into k subsets for k-fold cross validation.

# Arguments
- `spinconfig_list::AbstractVector{T}`: The dataset to be split.
- `k_fold::Integer`: The number of folds (subsets).
- `shuffle_data::Bool`: Whether to shuffle the data before splitting (default: true).

# Returns
- `Vector{Vector{T}}`: A vector of k subsets, each being a vector of elements from the original dataset.

# Throws
- `ArgumentError`: If `k_fold` is less than 2 or greater than the dataset size.
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
