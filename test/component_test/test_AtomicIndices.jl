using .CountingContainers
using .AtomicIndices

@testset "AtomicIndices" begin
	@testset "Indices" begin
		@testset "constructor" begin
			# Normal cases
			@test Indices(1, 1, 0, 1) isa Indices
			@test Indices(1, 2, -2, 1) isa Indices
			@test Indices(1, 2, 2, 1) isa Indices
		end

		@testset "domain errors" begin
			# l value errors
			@test_throws DomainError Indices(1, 0, -1, 1)
			@test_throws DomainError Indices(1, -1, 0, 1)
			
			# m value errors
			@test_throws DomainError Indices(1, 1, 2, 1)
			@test_throws DomainError Indices(1, 1, -2, 1)
			
			# cell value errors
			@test_throws DomainError Indices(1, 1, 0, 30)
			@test_throws DomainError Indices(1, 1, 0, 0)
			@test_throws DomainError Indices(1, 1, 0, -1)
		end

		@testset "comparison" begin
			# Equality
			@test Indices(1, 1, 0, 1) == Indices(1, 1, 0, 1)
			
			# Atom comparison
			@test Indices(1, 1, 0, 1) < Indices(2, 1, 0, 1)
			
			# l value comparison
			@test Indices(1, 1, 0, 1) < Indices(1, 2, 0, 1)
			
			# m value comparison
			@test Indices(1, 1, -1, 1) < Indices(1, 1, 0, 1)
			@test Indices(1, 1, 0, 1) < Indices(1, 1, 1, 1)
			
			# Cell comparison
			@test Indices(1, 1, 0, 1) < Indices(1, 1, 0, 2)
		end

		@testset "hash" begin
			# Same values should have same hash
			@test hash(Indices(1, 1, 0, 1)) == hash(Indices(1, 1, 0, 1))
			# Different values should have different hash
			@test hash(Indices(1, 1, 0, 1)) != hash(Indices(1, 1, 1, 1))
		end
	end

	@testset "IndicesUniqueList" begin
		@testset "constructor" begin
			# Empty list
			iul = IndicesUniqueList()
			@test isempty(iul)
			@test length(iul) == 0
			
			# Single element list
			indices = Indices(1, 1, 0, 1)
			iul = IndicesUniqueList(indices)
			@test length(iul) == 1
			@test indices in iul
			
			# Multiple elements list
			indices1 = Indices(1, 1, 0, 1)
			indices2 = Indices(2, 1, 0, 1)
			iul = IndicesUniqueList([indices1, indices2])
			@test length(iul) == 2
			@test indices1 in iul
			@test indices2 in iul
		end

		@testset "operations" begin
			indices1 = Indices(1, 1, 0, 1)
			indices2 = Indices(2, 1, 0, 1)
			indices3 = Indices(3, 1, 0, 1)
			
			iul = IndicesUniqueList()
			
			# push! operation
			push!(iul, indices1)
			@test length(iul) == 1
			@test indices1 in iul
			
			# Duplicate element addition
			push!(iul, indices1)
			@test length(iul) == 1
			
			# append! operation
			append!(iul, [indices2, indices3])
			@test length(iul) == 3
			@test indices2 in iul
			@test indices3 in iul
			
			# append! with duplicate elements
			append!(iul, [indices1, indices2])
			@test length(iul) == 3
		end

		@testset "getters" begin
			indices1 = Indices(1, 1, 0, 1)
			indices2 = Indices(2, 2, 0, 1)
			iul = IndicesUniqueList([indices1, indices2])
			
			@test get_total_L(iul) == 3
			@test get_atom_l_list(iul) == [[1, 1], [2, 2]]
		end

		@testset "product_indices_of_all_comb" begin
			# Simple case
			result = product_indices_of_all_comb([1], [1], [1])
			@test length(result) == 3  # m = -1, 0, 1
			
			# Multiple atoms case
			result = product_indices_of_all_comb([1, 2], [1, 1], [1, 1])
			@test length(result) == 9  # 3 * 3
			
			# Different l values case
			result = product_indices_of_all_comb([1, 2], [1, 2], [1, 1])
			@test length(result) == 24  # 3 * (3+5)

			# Error case: different length inputs
			@test_throws ErrorException product_indices_of_all_comb([1, 2], [1], [1, 1])
		end

		@testset "product_indices" begin
			# Simple case
			result = product_indices([1], [1], [1])
			@test length(result) == 3  # m = -1, 0, 1
			
			# Multiple atoms case
			result = product_indices([1, 2], [1, 1], [1, 1])
			@test length(result) == 9  # 3 * 3
			
			# Different l values case
			result = product_indices([1, 2], [1, 2], [1, 1])
			@test length(result) == 15  # 3 * 5

			# Error case: different length inputs
			@test_throws ErrorException product_indices([1, 2], [1], [1, 1])

			# Verify the order of combinations
			result = product_indices([1, 2], [1, 1], [1, 1])
			expected = [
				IndicesUniqueList([Indices(1, 1, -1, 1), Indices(2, 1, -1, 1)]),
				IndicesUniqueList([Indices(1, 1, -1, 1), Indices(2, 1, 0, 1)]),
				IndicesUniqueList([Indices(1, 1, -1, 1), Indices(2, 1, 1, 1)]),
				IndicesUniqueList([Indices(1, 1, 0, 1), Indices(2, 1, -1, 1)]),
				IndicesUniqueList([Indices(1, 1, 0, 1), Indices(2, 1, 0, 1)]),
				IndicesUniqueList([Indices(1, 1, 0, 1), Indices(2, 1, 1, 1)]),
				IndicesUniqueList([Indices(1, 1, 1, 1), Indices(2, 1, -1, 1)]),
				IndicesUniqueList([Indices(1, 1, 1, 1), Indices(2, 1, 0, 1)]),
				IndicesUniqueList([Indices(1, 1, 1, 1), Indices(2, 1, 1, 1)])
			]
			@test result == expected
		end

		@testset "indices_singleatom" begin
			# l=1 case
			result = indices_singleatom(1, 1, 1)
			@test length(result) == 3
			@test result == [
				Indices(1, 1, -1, 1),
				Indices(1, 1, 0, 1),
				Indices(1, 1, 1, 1)
			]
			
			# l=2 case
			result = indices_singleatom(1, 2, 1)
			@test length(result) == 5
			@test result == [
				Indices(1, 2, -2, 1),
				Indices(1, 2, -1, 1),
				Indices(1, 2, 0, 1),
				Indices(1, 2, 1, 1),
				Indices(1, 2, 2, 1)
			]
		end
	end
end

