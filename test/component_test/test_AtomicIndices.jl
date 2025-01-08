using .CountingContainers
using .AtomicIndices

@testset "AtomicIndices" begin
	@testset "Indices" begin
		indices1 = Indices(3, 1, -1)
		indices2 = Indices(4, 2, 1)
		indices3 = Indices(4, 2, -1)
		indices4 = Indices(4, 1, 0)
		@test indices1 == indices1
		@test indices1 == (3, 1, -1)
		@test (3, 1, -1) == indices1
		@test indices1 != indices2
		@test indices3 < indices2
		@test indices1 < indices2
		@test indices4 < indices3
	end

	@testset "IndicesUniqueList" begin
		indices1 = Indices(3, 1, -1)
		indices2 = Indices(4, 2, 1)
		indices3 = Indices(4, 2, -1)

		iul1 = IndicesUniqueList([indices1, indices2])
		@test length(iul1) == 2
		@test eltype(iul1) == Indices
		for val in iul1
		end
		@test indices1 in iul1
		@test !(indices3 in iul1)
		@test get_atomlist(iul1) == [3, 4]
		@test get_atomlist(iul1) != [4, 4]
		@test get_totalL(iul1) == 3
		@test get_totalL(iul1) != 4
		@test iul1 == iul1
		iul2 = IndicesUniqueList([indices1, indices3])
		@test iul2 < iul1

		iul3 = IndicesUniqueList()
		@test isempty(iul3)
		push!(iul3, indices2)
		push!(iul3, indices2)# ignored this operation because of the uniqueness
		@test iul3 < iul1
		push!(iul3, indices1)
		@test iul1 < iul3
		@test !(isempty(iul3))
		iul3_copy = copy(iul3)
		@test iul3 == iul3_copy
		@test iul3 != iul1
		@test sort(iul3) == iul1
		@test iul3[1] == indices2
		@test iul3[2] == indices1

		append!(iul3, iul2)
		@test iul3 == IndicesUniqueList([indices2, indices1, indices3])
		@test length(iul3) == 3

		@test indices_singleatom(3, 1) ==
			  [Indices(3, 1, -1), Indices(3, 1, 0), Indices(3, 1, 1)]
		@test indices_singleatom(4, 2) ==
			  [
			Indices(4, 2, -2),
			Indices(4, 2, -1),
			Indices(4, 2, 0),
			Indices(4, 2, 1),
			Indices(4, 2, 2),
		]

		iul_product = product_indices([4, 5], [1, 1])
		@test length(iul_product) == 3 * 3
		@test iul_product == [
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 1)]),
		]

		iul_product = product_indices([4, 5], [1, 2])
		@test length(iul_product) == 3 * (3 + 5)
		@test iul_product == [
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 2, -2)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 2, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 2, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 2, 1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 2, 2)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 2, -2)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 2, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 2, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 2, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 2, 2)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 2, -2)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 2, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 2, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 2, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 2, 2)]),
		]

		iul_product = product_indices([4, 5, 6], [1, 1, 1])
		@test length(iul_product) == 3 * 3 * 3
		@test iul_product == [
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, -1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, -1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, -1), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 0), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 0), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 0), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, -1), Indices(5, 1, 1), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, -1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, -1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, -1), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 0), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 0), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 0), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 0), Indices(5, 1, 1), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, -1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, -1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, -1), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 0), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 0), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 0), Indices(6, 1, 1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 1), Indices(6, 1, -1)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 1), Indices(6, 1, 0)]),
			IndicesUniqueList([Indices(4, 1, 1), Indices(5, 1, 1), Indices(6, 1, 1)])]

	end

end

