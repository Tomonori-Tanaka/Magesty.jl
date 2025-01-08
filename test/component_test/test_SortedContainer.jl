
using .SortedContainer

@testset "SoterdContainer" begin
	# ─────────────────────────────────────────────────────────────────────────────
	# SortedCountingUniqueVector
	# ─────────────────────────────────────────────────────────────────────────────
	vec = Vector{Int}([3, 3, 2, 1, 1])
	scv = SortedCountingUniqueVector(vec)
	@test scv == [1, 2, 3]
	@test getcount(scv, 1) == 2
	@test getcount(scv, 2) == 1
	@test getcount(scv, 3) == 2
	@test_throws KeyError scv.counts[0]


	# test push!
	push!(scv, 1)
	@test scv == [1, 2, 3]
	@test getcount(scv, 1) == 3

	# test append!
	append!(scv, [2, 3, 4])
	@test scv == [1, 2, 3, 4]
	@test getcount(scv, 1) == 3
	@test getcount(scv, 4) == 1
	@test getcount(scv, 2) == 2

	# test delete
	delete!(scv, 4)
	@test_throws KeyError scv.counts[4]
	@test scv == [1, 2, 3]

	# test in
	@test 1 in scv
	@test (5 in scv) == false

	# test copy and isless
	scv_copy = copy(scv)
	push!(scv_copy, 1)
	@test getcount(scv_copy, 1) != getcount(scv, 1)
	@test scv_copy.data == scv.data
	@test scv_copy != scv
	@test !(scv < scv_copy)


	# ─────────────────────────────────────────────────────────────────────────────
	# CountingUniqueVector
	# ─────────────────────────────────────────────────────────────────────────────

end
