using .AngularMomentumCoupling
using ..SphericalHarmonicsTransforms
using LinearAlgebra
using WignerSymbols
using OffsetArrays


@testset "AngularMomentumCoupling" begin
	@testset "mrange" begin
		@test mrange(0) == [0]
		@test mrange(1) == [-1, 0, 1]
		@test mrange(2) == [-2, -1, 0, 1, 2]
	end

	@testset "enumerate_paths_left_all" begin
		@test enumerate_paths_left_all([1]) == [(Int[], 1)]
		# l1 = l2 = 1
		@test enumerate_paths_left_all([1, 1]) == [(Int[], 0), (Int[], 1), (Int[], 2)]
		# l1 = 1, l2 = 2
		@test enumerate_paths_left_all([1, 2]) == [(Int[], 1), (Int[], 2), (Int[], 3)]
		# l1 = 1, l2 = 1, l3 = 1
		@test enumerate_paths_left_all([1, 1, 1]) ==
			  [
			(Int[0], 1),
			(Int[1], 0),
			(Int[1], 1),
			(Int[1], 2),
			(Int[2], 1),
			(Int[2], 2),
			(Int[2], 3),
		]
		# l1 = 1, l2 = 2, l3 = 1
		@test enumerate_paths_left_all([1, 2, 1]) ==
			  [
			(Int[1], 0),
			(Int[1], 1),
			(Int[1], 2),
			(Int[2], 1),
			(Int[2], 2),
			(Int[2], 3),
			(Int[3], 2),
			(Int[3], 3),
			(Int[3], 4),
		]
		# l1 = 1, l2 = 1, l3 = 1, l4 = 1
		@test enumerate_paths_left_all([1, 1, 1, 1]) ==
			  [
			(Int[0, 1], 0),
			(Int[0, 1], 1),
			(Int[0, 1], 2),
			(Int[1, 0], 1),
			(Int[1, 1], 0),
			(Int[1, 1], 1),
			(Int[1, 1], 2),
			(Int[1, 2], 1),
			(Int[1, 2], 2),
			(Int[1, 2], 3),
			(Int[2, 1], 0),
			(Int[2, 1], 1),
			(Int[2, 1], 2),
			(Int[2, 2], 1),
			(Int[2, 2], 2),
			(Int[2, 2], 3),
			(Int[2, 3], 2),
			(Int[2, 3], 3),
			(Int[2, 3], 4),
		]
	end
	@testset "coeff_tensor_complex" begin
		# l1 = 1, l2 = 1
		ls = [1, 1]
		paths = enumerate_paths_left_all(ls)
		paths_1 = paths[1] # (Lseq = [], Lf = 0)
		paths_2 = paths[2] # (Lseq = [], Lf = 1)
		paths_3 = paths[3] # (Lseq = [], Lf = 2)

		C = coeff_tensor_complex(ls, paths_1[1], paths_1[2])
		@test size(C) == (3, 3, 1)  # 1 = 2*Lf + 1, where Lf = 0
		# Lf = 0, Mf = 0: test as a 3x3 matrix
		Lf = paths_1[2]
		Mf = 0
		m1_range = mrange(ls[1])  # [-1, 0, 1]
		m2_range = mrange(ls[2])  # [-1, 0, 1]
		expected = zeros(Float64, 3, 3)
		for (i, m1) in enumerate(m1_range)
			for (j, m2) in enumerate(m2_range)
				expected[i, j] = Float64(clebschgordan(ls[1], m1, ls[2], m2, Lf, Mf))
			end
		end
		@test C[:, :, 1] ≈ expected

		C = coeff_tensor_complex(ls, paths_2[1], paths_2[2])
		@test size(C) == (3, 3, 3)  # 3 = 2*Lf + 1, where Lf = 1
		# Lf = 1: test each Mf slice (-1, 0, 1) as a 3x3 matrix
		Lf = paths_2[2]
		m1_range = mrange(ls[1])  # [-1, 0, 1]
		m2_range = mrange(ls[2])  # [-1, 0, 1]
		Mf_range = mrange(Lf)      # [-1, 0, 1]
		for (q, Mf) in enumerate(Mf_range)
			expected = zeros(Float64, 3, 3)
			for (i, m1) in enumerate(m1_range)
				for (j, m2) in enumerate(m2_range)
					expected[i, j] = Float64(clebschgordan(ls[1], m1, ls[2], m2, Lf, Mf))
				end
			end
			@test C[:, :, q] ≈ expected
		end

		C = coeff_tensor_complex(ls, paths_3[1], paths_3[2])
		@test size(C) == (3, 3, 5)  # 3 = 2*Lf + 1, where Lf = 2
		# Lf = 2: test each Mf slice (-2, -1, 0, 1, 2) as a 3x3 matrix
		Lf = paths_3[2]
		Mf_range = mrange(Lf)      # [-2, -1, 0, 1, 2]
		for (q, Mf) in enumerate(Mf_range)
			expected = zeros(Float64, 3, 3)
			for (i, m1) in enumerate(m1_range)
				for (j, m2) in enumerate(m2_range)
					expected[i, j] = Float64(clebschgordan(ls[1], m1, ls[2], m2, Lf, Mf))
				end
			end
			@test C[:, :, q] ≈ expected
		end
	end

	ls = [1, 1, 1]
	paths = enumerate_paths_left_all(ls)
	@test paths == [
		(Int[0], 1),
		(Int[1], 0),
		(Int[1], 1),
		(Int[1], 2),
		(Int[2], 1),
		(Int[2], 2),
		(Int[2], 3),
	]

	for (Lseq, Lf) in paths
		C = coeff_tensor_complex(ls, Lseq, Lf)
		@test size(C) == (3, 3, 3, 2*Lf + 1)
		m1_range = mrange(ls[1])
		m2_range = mrange(ls[2])
		m3_range = mrange(ls[3])
		Mf_range = mrange(Lf)
		expected = zeros(Float64, 3, 3, 3, 2*Lf + 1)
		for (i, m1) in enumerate(m1_range)
			for (j, m2) in enumerate(m2_range)
				for (k, m3) in enumerate(m3_range)
					for (q, Mf) in enumerate(Mf_range)
						if abs(m1+m2) > Lseq[1]
							continue
						end
						expected[i, j, k, q] =
							Float64(clebschgordan(ls[1], m1, ls[2], m2, Lseq[1], m1+m2)) *
							Float64(clebschgordan(Lseq[1], m1+m2, ls[3], m3, Lf, Mf))
					end
				end
			end
		end
		@test C ≈ expected
	end

	ls = [2, 2, 2]
	paths = enumerate_paths_left_all(ls)

	for (Lseq, Lf) in paths
		C = coeff_tensor_complex(ls, Lseq, Lf)
		@test size(C) == (5, 5, 5, 2*Lf + 1)
		m1_range = mrange(ls[1])
		m2_range = mrange(ls[2])
		m3_range = mrange(ls[3])
		Mf_range = mrange(Lf)
		expected = zeros(Float64, 5, 5, 5, 2*Lf + 1)
		for (i, m1) in enumerate(m1_range)
			for (j, m2) in enumerate(m2_range)
				for (k, m3) in enumerate(m3_range)
					for (q, Mf) in enumerate(Mf_range)
						if abs(m1+m2) > Lseq[1]
							continue
						end
						expected[i, j, k, q] =
							Float64(clebschgordan(ls[1], m1, ls[2], m2, Lseq[1], m1+m2)) *
							Float64(clebschgordan(Lseq[1], m1+m2, ls[3], m3, Lf, Mf))
					end
				end
			end
		end
		@test C ≈ expected
	end

	ls = [2, 2, 2, 2]
	paths = enumerate_paths_left_all(ls)
	for (Lseq, Lf) in paths
		C = coeff_tensor_complex(ls, Lseq, Lf)
		@test size(C) == (5, 5, 5, 5, 2*Lf + 1)
		m1_range = mrange(ls[1])
		m2_range = mrange(ls[2])
		m3_range = mrange(ls[3])
		m4_range = mrange(ls[4])
		Mf_range = mrange(Lf)
		expected = zeros(Float64, 5, 5, 5, 5, 2*Lf + 1)
		for (i, m1) in enumerate(m1_range)
			for (j, m2) in enumerate(m2_range)
				for (k, m3) in enumerate(m3_range)
					for (l, m4) in enumerate(m4_range)
						for (q, Mf) in enumerate(Mf_range)
							if abs(m1+m2) > Lseq[1]
								continue
							end
							if abs(m1+m2+m3) > Lseq[2]
								continue
							end
							expected[i, j, k, l, q] =
								Float64(clebschgordan(ls[1], m1, ls[2], m2, Lseq[1], m1+m2)) *
								Float64(
									clebschgordan(Lseq[1], m1+m2, ls[3], m3, Lseq[2], m1+m2+m3),
								) *
								Float64(clebschgordan(Lseq[2], m1+m2+m3, ls[4], m4, Lf, Mf))
						end
					end
				end
			end
		end
		@test C ≈ expected
	end

	@testset "build_all_real_bases" begin
		ls = [1, 1]
		bases, paths = build_all_real_bases(ls)

		L0_matrix = bases[0][1]
		L1_tensor = bases[1][1]
		L2_matrix = bases[2][1]
		L1_Mm1_matrix = L1_tensor[:, :, 1]
		L1_M0_matrix = L1_tensor[:, :, 2]
		L1_M1_matrix = L1_tensor[:, :, 3]
		L2_Mm2_matrix = L2_matrix[:, :, 1]
		L2_Mm1_matrix = L2_matrix[:, :, 2]
		L2_M0_matrix = L2_matrix[:, :, 3]
		L2_M1_matrix = L2_matrix[:, :, 4]
		L2_M2_matrix = L2_matrix[:, :, 5]

		# check normality
		@test sum((L0_matrix) .^ 2) ≈ 1.0
		@test sum((L1_Mm1_matrix) .^ 2) ≈ 1.0
		@test sum((L1_M0_matrix) .^ 2) ≈ 1.0
		@test sum((L1_M1_matrix) .^ 2) ≈ 1.0
		@test sum((L2_Mm2_matrix) .^ 2) ≈ 1.0
		@test sum((L2_Mm1_matrix) .^ 2) ≈ 1.0
		@test sum((L2_M0_matrix) .^ 2) ≈ 1.0
		@test sum((L2_M1_matrix) .^ 2) ≈ 1.0
		@test sum((L2_M2_matrix) .^ 2) ≈ 1.0

		# check orthogonality
		@test isapprox(sum(L0_matrix .* L1_Mm1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L1_M0_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L1_M1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L2_Mm2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L2_Mm1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L2_M0_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L2_M1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L0_matrix .* L2_M2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_Mm1_matrix .* L2_Mm1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_Mm1_matrix .* L2_M0_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_Mm1_matrix .* L2_M1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_Mm1_matrix .* L2_M2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M0_matrix .* L2_Mm2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M0_matrix .* L2_Mm1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M0_matrix .* L2_M0_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M0_matrix .* L2_M1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M0_matrix .* L2_M2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M1_matrix .* L2_Mm1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M1_matrix .* L2_M0_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M1_matrix .* L2_M1_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L1_M1_matrix .* L2_M2_matrix), 0.0, atol = 1e-10)
		@test isapprox(sum(L2_Mm2_matrix .* L2_Mm1_matrix), 0.0, atol = 1e-10)

		# display(bases[1][1][:, :, 3])

		# bases, paths = build_all_complex_bases(ls)
		# display(bases[1][1][:, :, 1])




		ls = [2, 1]
		bases, paths = build_all_complex_bases(ls)
		@test length(bases) == 3
		@test (length(bases[1]) == 1) && (length(bases[2]) == 1) && (length(bases[3]) == 1)
		@test length(paths) == 3
		@test (length(paths[1]) == 1) && (length(paths[2]) == 1) && (length(paths[3]) == 1)
		@test (paths[1][1] == Int[]) && (paths[2][1] == Int[]) && (paths[3][1] == Int[])

		# dimensionality check
		@test size(bases[1][1]) == (5, 3, 3)
		@test size(bases[2][1]) == (5, 3, 5)
		@test size(bases[3][1]) == (5, 3, 7)

		# Lf = 1
		Lf = 1
		tensor = bases[1][1]

		@test isapprox(tensor[:, :, 1],
			[ 0.0 0.0 sqrt(3/5);
				0.0 -sqrt(3/10) 0.0;
				sqrt(1/10) 0.0 0.0;
				0.0 0.0 0.0;
				0.0 0.0 0.0], atol = 1e-10)



	end
end