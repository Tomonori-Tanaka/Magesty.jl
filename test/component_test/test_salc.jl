using Test
using LinearAlgebra
using EzXML
using ..SALCs
using ..AtomicIndices

@testset "SALC Tests" begin
	@testset "Basic SALC Construction" begin
		# Test basic SALC construction
		basisset = Vector([
			IndicesUniqueList([Indices(1, 2, -2, 1), Indices(2, 3, 3, 4)]),
			IndicesUniqueList([Indices(1, 2, 1, 1), Indices(3, 2, 0, 1)]),
		])
		coeffs = [1.0, -1.0]
		coeffs = coeffs ./ norm(coeffs)
		multiplicity_list = [2, 1]

		salc = SALC(basisset, coeffs, multiplicity_list)

		@test length(salc.basisset) == 2
		@test length(salc.coeffs) == 2
		@test length(salc.multiplicity) == 2
		@test isapprox(salc.coeffs[1], 0.7071067811865475, atol = 1e-10)
		@test salc.multiplicity[1] == 2
		@test isapprox(salc.coeffs[2], -0.7071067811865475, atol = 1e-10)
		@test salc.multiplicity[2] == 1
	end

	@testset "SALC XML Initialization" begin
		@testset "Single Index Pattern" begin
			# XML string with single index pattern
			xml_str = """
			<SALC index="1">
				<basis multiplicity="2" index-1="1 2 -1 1">0.7071067811865475</basis>
				<basis multiplicity="1" index-1="3 4 2 1">-0.7071067811865475</basis>
			</SALC>
			"""

			# Parse XML document
			doc = parsexml(xml_str)
			salc = SALC(root(doc))

			# Test basic structure
			@test length(salc.basisset) == 2
			@test length(salc.coeffs) == 2
			@test length(salc.multiplicity) == 2

			# Test first basis vector
			@test salc.multiplicity[1] == 2
			@test isapprox(salc.coeffs[1], 0.7071067811865475, atol = 1e-10)  # 1.0/sqrt(2)
			@test length(salc.basisset[1]) == 1
			@test salc.basisset[1][1] == Indices(1, 2, -1, 1)

			# Test second basis vector
			@test salc.multiplicity[2] == 1
			@test isapprox(salc.coeffs[2], -0.7071067811865475, atol = 1e-10)  # -1.0/sqrt(2)
			@test length(salc.basisset[2]) == 1
			@test salc.basisset[2][1] == Indices(3, 4, 2, 1)
		end

		@testset "Double Index Pattern" begin
			# XML string with double index pattern
			xml_str = """
			<SALC index="1">
				<basis multiplicity="2" index-1="1 2 1 1" index-2="3 4 -4 6">0.7071067811865475</basis>
				<basis multiplicity="1" index-1="5 6 -1 1" index-2="7 8 0 10">-0.7071067811865475</basis>
			</SALC>
			"""

			# Parse XML document
			doc = parsexml(xml_str)
			salc = SALC(root(doc))

			# Test basic structure
			@test length(salc.basisset) == 2
			@test length(salc.coeffs) == 2
			@test length(salc.multiplicity) == 2

			# Test first basis vector
			@test salc.multiplicity[1] == 2
			@test isapprox(salc.coeffs[1], 0.7071067811865475, atol = 1e-10)  # 1.0/sqrt(2)
			@test length(salc.basisset[1]) == 2
			@test salc.basisset[1][1] == Indices(1, 2, 1, 1)
			@test salc.basisset[1][2] == Indices(3, 4, -4, 6)

			# Test second basis vector
			@test salc.multiplicity[2] == 1
			@test isapprox(salc.coeffs[2], -0.7071067811865475, atol = 1e-10)  # -1.0/sqrt(2)
			@test length(salc.basisset[2]) == 2
			@test salc.basisset[2][1] == Indices(5, 6, -1, 1)
			@test salc.basisset[2][2] == Indices(7, 8, 0, 10)
		end
	end
end
