module CrossValidation

using Random

include("../src/Magesty.jl")

export cross_validation
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
	spincluster::SpinCluster,
	weight_min::Real,
	weight_max::Real,
	weight_step::Integer,
	k_fold::Integer,
	;
	shuffle_data::Bool = true,
)
	# Parameter validation
	if weight_min >= weight_max
		throw(ArgumentError("weight_min ($weight_min) must be less than weight_max ($weight_max)"))
	end

	if weight_step <= 0
		throw(ArgumentError("weight_step ($weight_step) must be a positive integer"))
	end

	if k_fold < 2
		throw(ArgumentError("k_fold ($k_fold) must be at least 2"))
	end

	spinconfig_list::Vector{SpinConfig} = spincluster.optimize.spinconfig_list.spinconfigs
	if k_fold > length(spinconfig_list)
		throw(ArgumentError("k_fold ($k_fold) cannot be larger than the dataset size ($(length(spinconfig_list)))"))
	end

	splitted_spinconfig_list::Vector{Vector{SpinConfig}} =
		split_data(spinconfig_list, k_fold; shuffle_data = shuffle_data)
	test_rmse_list::Vector{Float64} = Vector{Float64}(undef, k_fold)
	train_rmse_list::Vector{Float64} = Vector{Float64}(undef, k_fold)
	
	current_spincluster::SpinCluster = spincluster
	for (i, weight::Float64) in enumerate(range(weight_min, weight_max, length = weight_step))
		cv_score_test::Float64, cv_score_train::Float64, current_spincluster::SpinCluster =
			cross_validation_with_a_weight(current_spincluster, weight, splitted_spinconfig_list)
		test_rmse_list[i] = cv_score_test
		train_rmse_list[i] = cv_score_train
	end

	optimum_weight::Float64 = range(weight_min, weight_max, length = weight_step)[argmin(test_rmse_list)]
	weight_list::Vector{Float64} = collect(range(weight_min, weight_max, length = weight_step))
	return optimum_weight, weight_list, test_rmse_list, train_rmse_list
end

function cross_validation_with_a_weight(
	spincluster::SpinCluster,
	weight::Real,
	splitted_dataset::Vector{Vector{SpinConfig}},
)
	test_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(splitted_dataset))
	train_rmse_list::Vector{Float64} = Vector{Float64}(undef, length(splitted_dataset))
	
	for (i, test_spinconfigs::Vector{SpinConfig}) in enumerate(splitted_dataset)
		# Organize training data
		train_spinconfigs::Vector{SpinConfig} =
			vcat(splitted_dataset[1:(i-1)]..., splitted_dataset[(i+1):end]...)

		# Create spin cluster for training data
		trained_spin_cluster::SpinCluster = SpinCluster(spincluster, weight, train_spinconfigs)
		# Calculate RMSE for training data
		predicted_energies_train::Vector{Float64} = trained_spin_cluster.optimize.predicted_energy_list
		observed_energies_train::Vector{Float64} =
			[train_spinconfig.energy for train_spinconfig in train_spinconfigs]
		train_rmse_list[i] =
			sqrt(mean((predicted_energies_train .- observed_energies_train) .^ 2))

		# Calculate RMSE for test data
		predicted_energies_test::Vector{Float64} = Vector{Float64}(undef, length(test_spinconfigs))
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
	return cv_score_test, cv_score_train, trained_spin_cluster
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
