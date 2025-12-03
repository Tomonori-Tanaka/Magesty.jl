@testset "Basis" begin
	using .AngularMomentumCoupling
	using Magesty.BasisSets: classify_tesseral_basislist_test
	using .Basis: LinearCombo, tesseral_linear_combos_from_tesseral_bases
	using .SortedContainer: SortedCountingUniqueVector
	using Test

	@testset "LinearCombo" begin
		@testset "Constructor" begin
			atoms = [10, 20]
			ls = [1, 1]
			lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
			@test length(lc_list) == 9
			# Lf = 0
			lc_lf0 = lc_list[1]
			@test lc_lf0.ls == (1, 1)
			@test lc_lf0.Lf == 0
			@test lc_lf0.Lseq == Int[]
			@test lc_lf0.atoms == [10, 20]
			@test lc_lf0.coeff_list == [1.0]
			@test isapprox(
				lc_lf0.coeff_tensor,
				[1/sqrt(3) 0.0 0.0; 0.0 1/sqrt(3) 0.0; 0.0 0.0 1/sqrt(3)],
				atol = 1e-10,
			)
		end

		@testset "dot" begin
			using LinearAlgebra

			@testset "Same LinearCombo" begin
				atoms = [1, 2]
				ls = [1, 1]
				lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
				lc1 = lc_list[1]
				# Same LinearCombo should have non-zero dot product (coeff_list dot itself)
				for (i, lc1) in enumerate(lc_list), (j, lc2) in enumerate(lc_list)
					if i == j
						@test dot(lc1, lc2) ≈ 1.0
					else
						@test dot(lc1, lc2) == 0.0
					end
				end
				lc_lf1 = lc_list[findfirst(lc -> lc.Lf == 1, lc_list)]
				lc_lf1_2 = deepcopy(lc_lf1)
				coeff_list = [1/sqrt(2), 1/sqrt(2), 0.0]
				lc_lf1_2 = LinearCombo(
					lc_lf1.ls,
					lc_lf1.Lf,
					lc_lf1.Lseq,
					lc_lf1.atoms,
					coeff_list,
					lc_lf1.coeff_tensor,
				)
				@test dot(lc_lf1, lc_lf1_2) ≈ 1/sqrt(2)

			end

			@testset "Different Lf" begin
				atoms = [1, 2]
				ls = [1, 1]
				lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
				# Find LinearCombos with different Lf
				lc_lf0 = lc_list[findfirst(lc -> lc.Lf == 0, lc_list)]
				lc_lf1 = lc_list[findfirst(lc -> lc.Lf == 1, lc_list)]
				@test dot(lc_lf0, lc_lf1) == 0.0
			end

			@testset "Different Lseq" begin
				atoms = [1, 2, 3]
				ls = [1, 1, 1]
				lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
				# Find LinearCombos with same Lf but different Lseq
				lcs_lf1 = filter(lc -> lc.Lf == 1, lc_list)
				if length(lcs_lf1) >= 2
					lc1 = lcs_lf1[1]
					lc2 = lcs_lf1[2]
					if lc1.Lseq != lc2.Lseq
						@test dot(lc1, lc2) == 0.0
					end
				end
			end

			@testset "Permuted atoms (same ls)" begin
				atoms1 = [1, 2]
				atoms2 = [2, 1]
				ls = [1, 1]
				lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls, atoms1)
				lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls, atoms2)
				for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
					if i == j
						@test dot(lc1, lc2) ≈ 1.0
					else
						@test dot(lc1, lc2) == 0.0
					end
				end
			end

			@testset "Permuted atoms (different ls order but same multiset)" begin
				atoms1 = [1, 2]
				atoms2 = [2, 1]
				ls1 = [2, 2]
				ls2 = [2, 2]
				lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls1, atoms1)
				lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls2, atoms2)
				for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
					if i == j
						@test dot(lc1, lc2) ≈ 1.0
					else
						@test dot(lc1, lc2) == 0.0
					end
				end
			end

			@testset "Non-matching (atom, l) pairs" begin
				atoms1 = [1, 2]
				atoms2 = [1, 3]  # Different atom
				ls = [1, 1]
				Lf = 0
				Lseq = Int[]
				coeff_list = [1.0]
				coeff_tensor = [1/sqrt(3) 0.0 0.0; 0.0 1/sqrt(3) 0.0; 0.0 0.0 1/sqrt(3)]
				lc1 = LinearCombo(ls, Lf, Lseq, atoms1, coeff_list, coeff_tensor)
				lc2 = LinearCombo(ls, Lf, Lseq, atoms2, coeff_list, coeff_tensor)
				# Should be zero because (atom, l) pairs don't match
				@test dot(lc1, lc2) == 0.0
			end

			@testset "Three-body case" begin
				@testset "ls=[1,1,1], atoms permuted" begin
					atoms1 = [1, 2, 3]
					atoms2 = [3, 1, 2]  # Permuted
					ls = [1, 1, 1]
					lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls, atoms1)
					lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls, atoms2)
					# Find matching Lf and Lseq
					lc1 = lc_list1[1]  # Lf=0
					lc2 = lc_list2[1]  # Lf=0
					# Should match because (atom, l) pairs are the same
					result = dot(lc1, lc2)
					for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
						if i == j
							@test dot(lc1, lc2) ≈ 1.0
						else
							@test dot(lc1, lc2) == 0.0
						end
					end
				end

				@testset "ls=[2,1,1], atoms=[3,1,2] vs [3,2,1]" begin
					atoms1 = [3, 1, 2]
					atoms2 = [3, 2, 1]
					ls = [2, 1, 1]
					lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls, atoms1)
					lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls, atoms2)
					# (atom, l) pairs: [(3,2), (1,1), (2,1)] vs [(3,2), (2,1), (1,1)]
					# These should match via permutation
					for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
						if lc1.Lf == lc2.Lf && lc1.Lseq == lc2.Lseq
							# Same Lf and Lseq: should have dot product equal to coeff_list dot product
							result = dot(lc1, lc2)
							expected = dot(lc1.coeff_list, lc2.coeff_list)
							@test result ≈ expected
						else
							# Different Lf or Lseq: should be zero
							@test dot(lc1, lc2) == 0.0
						end
					end
				end
			end

			@testset "Different number of sites" begin
				atoms1 = [1, 2]
				atoms2 = [1, 2, 3]
				ls1 = [1, 1]
				ls2 = [1, 1, 1]
				lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls1, atoms1)
				lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls2, atoms2)
				for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
					@test dot(lc1, lc2) == 0.0
				end
			end

			@testset "N=1 (single atom) case" begin
				@testset "Same LinearCombo" begin
					atoms = [1]
					ls = [2]
					lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
					@test length(lc_list) == 5  # 2*Lf+1 = 2*2+1 = 5
					# Check that Lseq is empty for N=1
					for lc in lc_list
						@test lc.Lseq == Int[]
						@test length(lc.coeff_list) == 5
						@test size(lc.coeff_tensor) == (5,)
					end
					# Same LinearCombo should have dot product of 1.0
					for (i, lc1) in enumerate(lc_list), (j, lc2) in enumerate(lc_list)
						if i == j
							@test dot(lc1, lc2) ≈ 1.0
						else
							@test dot(lc1, lc2) == 0.0
						end
					end
				end

				@testset "Different Lf" begin
					atoms = [1]
					ls = [2]
					lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
					# Find LinearCombos with different Lf
					lc_lf2 = lc_list[findfirst(lc -> lc.Lf == 2, lc_list)]
					# Create another LinearCombo with Lf=0
					ls_lf0 = [2]
					Lf0 = 0
					Lseq0 = Int[]
					coeff_list0 = [1.0]
					coeff_tensor0 = [1.0]
					lc_lf0 = LinearCombo(ls_lf0, Lf0, Lseq0, atoms, coeff_list0, coeff_tensor0)
					@test dot(lc_lf2, lc_lf0) == 0.0
				end

				@testset "Different atom" begin
					atoms1 = [1]
					atoms2 = [2]
					ls = [2]
					lc_list1 = tesseral_linear_combos_from_tesseral_bases(ls, atoms1)
					lc_list2 = tesseral_linear_combos_from_tesseral_bases(ls, atoms2)
					# Different atoms should give zero dot product
					for (i, lc1) in enumerate(lc_list1), (j, lc2) in enumerate(lc_list2)
						@test dot(lc1, lc2) == 0.0
					end
				end

				@testset "Different coeff_list" begin
					atoms = [1]
					ls = [2]
					lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
					lc1 = lc_list[1]
					# Create a LinearCombo with different coeff_list
					coeff_list2 = [1/sqrt(2), 1/sqrt(2), 0.0, 0.0, 0.0]
					lc2 = LinearCombo(
						lc1.ls,
						lc1.Lf,
						lc1.Lseq,
						lc1.atoms,
						coeff_list2,
						lc1.coeff_tensor,
					)
					# Dot product should be the dot product of coeff_list vectors
					expected = dot(lc1.coeff_list, coeff_list2)
					@test dot(lc1, lc2) ≈ expected
				end
			end
		end
	end

	@testset "classify_tesseral_basislist_test" begin
		atoms2 = [1, 2]
		atoms3 = [1, 2, 3]
		ls2 = [1, 1]
		ls3 = [1, 1, 1]
		two_body_list = tesseral_linear_combos_from_tesseral_bases(ls2, atoms2)
		three_body_list = tesseral_linear_combos_from_tesseral_bases(ls3, atoms3)

		lc_two_lf1 = filter(lc -> lc.Lf == 1, two_body_list)
		@test length(lc_two_lf1) >= 2
		lf0_idx_two = findfirst(lc -> lc.Lf == 0, two_body_list)
		@test !isnothing(lf0_idx_two)
		lf0_idx_three = findfirst(lc -> lc.Lf == 0, three_body_list)
		@test !isnothing(lf0_idx_three)

		lc_lf1_a = lc_two_lf1[1]
		lc_lf1_b = lc_two_lf1[2]
		lc_two_lf0 = two_body_list[lf0_idx_two]
		lc_three_lf0 = three_body_list[lf0_idx_three]

		sample_list = [lc_lf1_a, lc_lf1_b, lc_two_lf0, lc_three_lf0]
		result = classify_tesseral_basislist_test(sample_list)

		@test length(result) == 3
		@test length(result[1]) == 2
		@test all(lc -> length(lc.atoms) == 2 && lc.Lf == 1, result[1])
		@test length(result[2]) == 1
		@test only(result[2]).Lf == 0
		@test length(result[3]) == 1
		@test length(only(result[3]).atoms) == 3

		scv = SortedCountingUniqueVector{Basis.LinearCombo}()
		push!(scv, lc_lf1_a)
		push!(scv, lc_lf1_a)
		push!(scv, lc_three_lf0)
		counted_result = classify_tesseral_basislist_test(scv)

		@test counted_result[1].counts[lc_lf1_a] == 2
		@test counted_result[2].counts[lc_three_lf0] == 1
	end
end
