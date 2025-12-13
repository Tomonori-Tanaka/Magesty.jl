using .CountingContainer
using .AtomicIndices

@testset "AtomicIndices" begin
	@testset "shsiteindex_singleatom" begin
		result = shsiteindex_singleatom(1, 1)
		@test length(result) == 3
		@test result == [
			SHSiteIndex(1, 1, -1),
			SHSiteIndex(1, 1, 0),
			SHSiteIndex(1, 1, 1),
		]

		result = shsiteindex_singleatom(1, 2)
		@test length(result) == 5
		@test result == [
			SHSiteIndex(1, 2, -2),
			SHSiteIndex(1, 2, -1),
			SHSiteIndex(1, 2, 0),
			SHSiteIndex(1, 2, 1),
			SHSiteIndex(1, 2, 2),
		]
	end

	@testset "product_shsiteindex" begin
		result = product_shsiteindex([1], [1])
		@test length(result) == 3
		@test result == [
			SHProduct([SHSiteIndex(1, 1, -1)]),
			SHProduct([SHSiteIndex(1, 1, 0)]),
			SHProduct([SHSiteIndex(1, 1, 1)]),
		]

		result = product_shsiteindex([1, 2], [1, 1])
		@test length(result) == 9

		result = product_shsiteindex([1, 2], [1, 2])
		@test length(result) == 15
		expected = [
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 2, -2)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 2, -1)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 2, 0)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 2, 1)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 2, 2)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 2, -2)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 2, -1)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 2, 0)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 2, 1)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 2, -2)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 2, -1)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 2, 0)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 2, 1)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 2, 2)]),
		]

		result = product_shsiteindex([1, 2], [1, 1])
		expected = [
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 1, -1)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 1, 0)]),
			SHProduct([SHSiteIndex(1, 1, -1), SHSiteIndex(2, 1, 1)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 1, -1)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 1, 0)]),
			SHProduct([SHSiteIndex(1, 1, 0), SHSiteIndex(2, 1, 1)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 1, -1)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 1, 0)]),
			SHProduct([SHSiteIndex(1, 1, 1), SHSiteIndex(2, 1, 1)]),
		]
		@test result == expected

		result = product_shsiteindex([1, 2], [1, 2])
		@test length(result) == 15
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
			
			
			# append! operation
			append!(shp, [shsi2, shsi3])
			@test length(shp) == 3
			@test shsi2 in shp
			@test shsi3 in shp
			
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
	end
end

