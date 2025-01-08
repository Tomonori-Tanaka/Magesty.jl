using .SortedContainer
using .AtomicIndices

@testset "AtomicIndices" begin
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

	iul1 = IndicesSortedUniqueList([indices1, indices2])
	@test getatoms(iul1) == [3, 4]
	@test getatoms(iul1) != [4, 4]
	@test gettotall(iul1) == 3
	@test gettotall(iul1) != 4
	@test iul1 == iul1
	iul2 = IndicesSortedUniqueList([indices1, indices3])
	@test iul2 < iul1

	scv = SortedCountingUniqueVector([iul1, iul1, iul2])
	scv2 = SortedCountingUniqueVector([iul1, iul2, iul2])
	@test scv[1] == iul2
	@test scv[2] == iul1
	@test length(scv) == 2
	@test scv == sort([iul1, iul2])
	@test sort([iul1, iul2]) == scv
	@test scv == scv
	@test scv != scv2
	@test getcount(scv, iul1) == 2
	@test getcount(scv, iul2) == 1
	push!(scv, iul1)
	@test getcount(scv, iul1) == 3
	append!(scv2, [iul1, iul2])
	@test getcount(scv2, iul1) == 2
	@test getcount(scv2, iul2) == 3
	delete!(scv2, iul2)
	@test iul1 in scv2
	@test !(iul2 in scv2)

	scv3 = SortedCountingUniqueVector([iul1, iul1, iul2])
	scv3_copy = copy(scv3)
	push!(scv3_copy, iul1)
	@test scv3 != scv3_copy
end

