using .CountingContainers

@testset "CountingContainer" begin
	# ─────────────────────────────────────────────────────────────────────────────
	# CountingUniqueVector
	# ───────────────────────────────────────────────────────────────────────────── 

    # Int type case
	vec = Vector{Int}([3, 3, 2, 1, 1])
	cuv1 = CountingUniqueVector(vec)
	cuv1_1 = CountingUniqueVector{Int}()
	for (idx, elm) in enumerate(vec)
		push!(cuv1_1, elm)
	end
	@test cuv1 == cuv1_1
	cuv1_2 = CountingUniqueVector{Int}()
	append!(cuv1_2, vec)
	@test cuv1 == cuv1_2
	@test length(cuv1) == length(cuv1_1)
	for (idx, elm) in enumerate(cuv1)
		@test elm == cuv1_1[idx]
		@test elm == cuv1_2[idx]
	end

	@test getcounts(cuv1, 3) == 2
	@test getcounts(cuv1, 2) == 1
	@test getcounts(cuv1, 1) == 2
	@test getcounts(cuv1, 4) == 0
	@test size(cuv1) == (3,)
	@test 3 in cuv1
	@test 2 in cuv1
	@test 1 in cuv1
	@test !(4 in cuv1)
	@test_throws KeyError cuv1.counts[4]
	cuv1_copy = copy(cuv1)
	@test cuv1_copy == cuv1


end
