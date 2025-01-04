
using .AtomicIndices

@testset "AtomicIndices" begin
	indices1 = Indices(3, 2, 1, -1)
	indices2 = Indices(4, 3, 2, 1)
	indices3 = Indices(4, 3, 2, -1)
	indices4 = Indices(4, 3, 1, 0)
	@test indices1 == indices1
	@test indices1 == (3, 2, 1, -1)
	@test (3, 2, 1, -1) == indices1
	@test indices1 != indices2
	@test indices3 < indices2
	@test indices1 < indices2
	@test indices4 < indices3

	iul1 = IndicesUniqueList([indices1, indices2])
	@test getatoms(iul1) == (3, 4)
	@test getatoms(iul1) != (4, 4)
	@test gettotall(iul1) == 3
	@test gettotall(iul1) != 4
	@test iul1 == iul1
	iul2 = IndicesUniqueList([indices1, indices3])
	@test iul2 < iul1
end
