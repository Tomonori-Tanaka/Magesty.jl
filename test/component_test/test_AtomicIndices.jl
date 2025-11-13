using .CountingContainer
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

	@testset "SHSiteIndex" begin
		@testset "constructor" begin
			# Normal cases
			@test SHSiteIndex(1, 1, 0) isa SHSiteIndex
			@test SHSiteIndex(1, 2, -2) isa SHSiteIndex
			@test SHSiteIndex(1, 2, 2) isa SHSiteIndex
			@test SHSiteIndex((1, 1, 0)) isa SHSiteIndex
		end

		@testset "domain errors" begin
			# l value errors
			@test_throws DomainError SHSiteIndex(1, 0, -1)
			@test_throws DomainError SHSiteIndex(1, -1, 0)
			
			# m value errors
			@test_throws DomainError SHSiteIndex(1, 1, 2)
			@test_throws DomainError SHSiteIndex(1, 1, -2)
		end

		@testset "comparison" begin
			# Equality
			@test SHSiteIndex(1, 1, 0) == SHSiteIndex(1, 1, 0)
			@test SHSiteIndex(1, 1, 0) == (1, 1, 0)
			@test (1, 1, 0) == SHSiteIndex(1, 1, 0)
			
			# i comparison
			@test SHSiteIndex(1, 1, 0) < SHSiteIndex(2, 1, 0)
			
			# l value comparison
			@test SHSiteIndex(1, 1, 0) < SHSiteIndex(1, 2, 0)
			
			# m value comparison
			@test SHSiteIndex(1, 1, -1) < SHSiteIndex(1, 1, 0)
			@test SHSiteIndex(1, 1, 0) < SHSiteIndex(1, 1, 1)
			
			# Tuple comparison
			@test SHSiteIndex(1, 1, 0) < (2, 1, 0)
			@test (1, 1, 0) < SHSiteIndex(2, 1, 0)
		end

		@testset "hash" begin
			# Same values should have same hash
			@test hash(SHSiteIndex(1, 1, 0)) == hash(SHSiteIndex(1, 1, 0))
			# Different values should have different hash
			@test hash(SHSiteIndex(1, 1, 0)) != hash(SHSiteIndex(1, 1, 1))
		end

		@testset "tuple conversion" begin
			# Convert to tuple
			shsi = SHSiteIndex(1, 2, -1)
			tuple_result = convert(NTuple{3, Int}, shsi)
			@test tuple_result == (1, 2, -1)
			
			# Convert from tuple
			shsi_from_tuple = SHSiteIndex((1, 2, -1))
			@test shsi_from_tuple == shsi
		end
	end

	@testset "SHProduct" begin
		@testset "constructor" begin
			# Empty product
			shp = SHProduct()
			@test isempty(shp)
			@test length(shp) == 0
			
			# Single element product
			shsi = SHSiteIndex(1, 1, 0)
			shp = SHProduct(shsi)
			@test length(shp) == 1
			@test shsi in shp
			
			# Multiple elements product
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			shp = SHProduct([shsi1, shsi2])
			@test length(shp) == 2
			@test shsi1 in shp
			@test shsi2 in shp
			
			# Duplicate elements are removed and sorted
			shsi3 = SHSiteIndex(3, 1, 0)
			shp = SHProduct([shsi3, shsi1, shsi2, shsi1])
			@test length(shp) == 3
			@test shp[1] == shsi1
			@test shp[2] == shsi2
			@test shp[3] == shsi3
		end

		@testset "operations" begin
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			shsi3 = SHSiteIndex(3, 1, 0)
			
			shp = SHProduct()
			
			# push! operation
			push!(shp, shsi1)
			@test length(shp) == 1
			@test shsi1 in shp
			
			# Duplicate element addition
			push!(shp, shsi1)
			@test length(shp) == 1
			
			# append! operation
			append!(shp, [shsi2, shsi3])
			@test length(shp) == 3
			@test shsi2 in shp
			@test shsi3 in shp
			
			# append! with duplicate elements
			append!(shp, [shsi1, shsi2])
			@test length(shp) == 3
		end

		@testset "comparison" begin
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			shsi3 = SHSiteIndex(3, 1, 0)
			
			# Equality
			shp1 = SHProduct([shsi1, shsi2])
			shp2 = SHProduct([shsi1, shsi2])
			@test shp1 == shp2
			
			# Different products
			shp3 = SHProduct([shsi1, shsi3])
			@test shp1 != shp3
			
			# Length comparison
			shp_single = SHProduct([shsi1])
			@test shp_single < shp1
			
			# Lexicographic comparison
			shp4 = SHProduct([shsi2, shsi3])
			@test shp1 < shp4
		end

		@testset "hash" begin
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			
			# Same products should have same hash
			shp1 = SHProduct([shsi1, shsi2])
			shp2 = SHProduct([shsi1, shsi2])
			@test hash(shp1) == hash(shp2)
			
			# Different products should have different hash
			shp3 = SHProduct([shsi1])
			@test hash(shp1) != hash(shp3)
		end

		@testset "AbstractVector interface" begin
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			shsi3 = SHSiteIndex(3, 1, 0)
			shp = SHProduct([shsi1, shsi2, shsi3])
			
			# getindex
			@test shp[1] == shsi1
			@test shp[2] == shsi2
			@test shp[3] == shsi3
			
			# length
			@test length(shp) == 3
			
			# eltype
			@test eltype(SHProduct) == SHSiteIndex
			
			# iterate
			collected = [x for x in shp]
			@test collected == [shsi1, shsi2, shsi3]
			
			# isempty
			@test !isempty(shp)
			@test isempty(SHProduct())
			
			# size
			@test size(shp) == (3,)
			
			# in
			@test shsi1 in shp
			@test shsi2 in shp
			@test !(SHSiteIndex(4, 1, 0) in shp)
		end

		@testset "sort" begin
			shsi1 = SHSiteIndex(1, 1, 0)
			shsi2 = SHSiteIndex(2, 1, 0)
			shsi3 = SHSiteIndex(3, 1, 0)
			
			# Create product using push! to control order (constructor auto-sorts)
			shp = SHProduct()
			push!(shp, shsi3)
			push!(shp, shsi1)
			push!(shp, shsi2)
			
			# Verify original order (unsorted)
			@test shp[1] == shsi3
			@test shp[2] == shsi1
			@test shp[3] == shsi2
			
			# Sort creates a new sorted product
			sorted_shp = sort(shp)
			
			@test sorted_shp[1] == shsi1
			@test sorted_shp[2] == shsi2
			@test sorted_shp[3] == shsi3
			
			# Original is not modified
			@test shp[1] == shsi3
			@test shp[2] == shsi1
			@test shp[3] == shsi2
		end
	end
end

