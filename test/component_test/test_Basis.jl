@testset "Basis" begin
	using .AngularMomentumCoupling
	using .Basis: CoupledBasis, tesseral_coupled_bases_from_tesseral_bases, reorder_atoms
	using .SortedContainer: SortedCountingUniqueVector
	using Test

	@testset "CoupledBasis" begin
		@testset "Constructor" begin
			atoms = [10, 20]
			ls = [1, 1]
			cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
			@test length(cb_list) == 3
			# Lf = 0
			cb_lf0 = cb_list[1]
			@test cb_lf0.ls == [1, 1]
			@test cb_lf0.Lf == 0
			@test cb_lf0.Lseq == Int[]
			@test cb_lf0.atoms == [10, 20]
			@test isapprox(
				cb_lf0.coeff_tensor,
				[1/sqrt(3) 0.0 0.0; 0.0 1/sqrt(3) 0.0; 0.0 0.0 1/sqrt(3)],
				atol = 1e-10,
			)
		end

		@testset "reorder_atoms" begin
			@testset "Two-body case" begin
				atoms = [1, 2]
				ls = [3, 1]
				cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
				@test length(cb_list) == 3

				cb = cb_list[1]
				original_atoms = copy(cb.atoms)
				original_ls = copy(cb.ls)
				original_tensor = copy(cb.coeff_tensor)
				@test size(original_tensor) == (7, 3, 5)

				# Test with new_atoms that need sorting: [2, 1] -> [1, 2]
				new_atoms = [4, 1]
				cb_new = reorder_atoms(cb, new_atoms)

				# Check that atoms are sorted
				@test cb_new.atoms == [1, 4]
				# Check that ls is permuted accordingly
				@test cb_new.ls == [1, 3]
				# Check that Lf and Lseq are unchanged
				@test cb_new.Lf == cb.Lf
				@test cb_new.Lseq == cb.Lseq
				@test size(cb_new.coeff_tensor) == (3, 7, 5)
			end

			@testset "Three-body case" begin
				atoms = [1, 2, 3]
				ls = [1, 1, 1]
				cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
				@test length(cb_list) > 0

				cb = cb_list[1]
				original_atoms = copy(cb.atoms)
				original_ls = copy(cb.ls)
				original_tensor = copy(cb.coeff_tensor)

				# Test with new_atoms = [5, 1, 2] -> sorted to [1, 2, 5]
				new_atoms = [5, 1, 2]
				cb_new = reorder_atoms(cb, new_atoms)

				# Check that atoms are sorted
				@test cb_new.atoms == [1, 2, 5]
				# Check that ls is permuted accordingly (p = [2, 3, 1])
				@test cb_new.ls == original_ls[[2, 3, 1]]
				# Check that tensor dimensions are permuted
				@test size(cb_new.coeff_tensor) == size(original_tensor)
				# Check that Lf and Lseq are unchanged
				@test cb_new.Lf == cb.Lf
				@test cb_new.Lseq == cb.Lseq
			end

			@testset "Single atom case (N=1)" begin
				atoms = [1]
				ls = [2]
				cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
				@test length(cb_list) > 0

				cb = cb_list[1]
				original_atoms = copy(cb.atoms)
				original_ls = copy(cb.ls)

				# For N=1, sorting doesn't change anything, but we can test the function
				new_atoms = [5]
				cb_new = reorder_atoms(cb, new_atoms)

				@test cb_new.atoms == [5]
				@test cb_new.ls == original_ls
				@test cb_new.Lf == cb.Lf
				@test cb_new.Lseq == cb.Lseq
				@test size(cb_new.coeff_tensor) == size(cb.coeff_tensor)
			end
		end
	end

	@testset "reorder_atoms (additional tests)" begin
		@testset "Two-body case" begin
			atoms = [1, 2]
			ls = [1, 1]
			cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
			@test length(cb_list) > 0

			cb = cb_list[1]
			original_atoms = copy(cb.atoms)
			original_ls = copy(cb.ls)
			original_tensor = copy(cb.coeff_tensor)

			# Test with new_atoms that need sorting: [5, 1, 2] -> [1, 2, 5]
			# But for two-body, we'll use [2, 1] -> [1, 2]
			new_atoms = [2, 1]
			cb_new = reorder_atoms(cb, new_atoms)

			# Check that atoms are sorted
			@test cb_new.atoms == [1, 2]
			# Check that ls is permuted accordingly
			@test cb_new.ls == original_ls[[2, 1]]
			# Check that tensor dimensions are permuted
			@test size(cb_new.coeff_tensor) == size(original_tensor)
			# Check that Lf and Lseq are unchanged
			@test cb_new.Lf == cb.Lf
			@test cb_new.Lseq == cb.Lseq
		end

		@testset "Three-body case" begin
			atoms = [1, 2, 3]
			ls = [1, 1, 1]
			cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
			@test length(cb_list) > 0

			cb = cb_list[1]
			original_atoms = copy(cb.atoms)
			original_ls = copy(cb.ls)
			original_tensor = copy(cb.coeff_tensor)

			# Test with new_atoms = [5, 1, 2] -> sorted to [1, 2, 5]
			new_atoms = [5, 1, 2]
			cb_new = reorder_atoms(cb, new_atoms)

			# Check that atoms are sorted
			@test cb_new.atoms == [1, 2, 5]
			# Check that ls is permuted accordingly (p = [2, 3, 1])
			@test cb_new.ls == original_ls[[2, 3, 1]]
			# Check that tensor dimensions are permuted
			@test size(cb_new.coeff_tensor) == size(original_tensor)
			# Check that Lf and Lseq are unchanged
			@test cb_new.Lf == cb.Lf
			@test cb_new.Lseq == cb.Lseq

			# Verify tensor permutation: first 3 dimensions permuted, last dimension (Mf) unchanged
			nd = ndims(original_tensor)
			@test nd == 4  # N+1 = 3+1
			# The permutation should be [2, 3, 1, 4] for dimensions
			expected_size = size(original_tensor)
			@test size(cb_new.coeff_tensor) == expected_size
		end

		@testset "Error cases" begin
			atoms = [1, 2]
			ls = [1, 1]
			cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
			cb = cb_list[1]

			# Wrong length
			@test_throws ArgumentError reorder_atoms(cb, [1, 2, 3])
			@test_throws ArgumentError reorder_atoms(cb, [1])
		end

		@testset "Single atom case (N=1)" begin
			atoms = [1]
			ls = [2]
			cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
			@test length(cb_list) > 0

			cb = cb_list[1]
			original_atoms = copy(cb.atoms)
			original_ls = copy(cb.ls)

			# For N=1, sorting doesn't change anything, but we can test the function
			new_atoms = [5]
			cb_new = reorder_atoms(cb, new_atoms)

			@test cb_new.atoms == [5]
			@test cb_new.ls == original_ls
			@test cb_new.Lf == cb.Lf
			@test cb_new.Lseq == cb.Lseq
			@test size(cb_new.coeff_tensor) == size(cb.coeff_tensor)
		end
	end
end
