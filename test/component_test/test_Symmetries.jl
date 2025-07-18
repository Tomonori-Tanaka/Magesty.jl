using Test
using StaticArrays
using Magesty.Symmetries
using Magesty.Structures

@testset "Symmetry Tests" begin
	@testset "Basic Components" begin
		@testset "SymmetryOperation Construction" begin
			# Test basic symmetry operation (identity)
			rotation_frac = SMatrix{3, 3, Float64}([
				1.0 0.0 0.0;
				0.0 1.0 0.0;
				0.0 0.0 1.0
			])
			rotation_cart = SMatrix{3, 3, Float64}([
				1.0 0.0 0.0;
				0.0 1.0 0.0;
				0.0 0.0 1.0
			])
			translation_frac = SVector{3, Float64}([0.0, 0.0, 0.0])

			symop = SymmetryOperation(
				rotation_frac,
				rotation_cart,
				translation_frac,
				false,  # is_translation
				false,  # is_translation_included
				true,   # is_proper
			)

			@test symop.rotation_frac == rotation_frac
			@test symop.rotation_cart == rotation_cart
			@test symop.translation_frac == translation_frac
			@test !symop.is_translation
			@test !symop.is_translation_included
			@test symop.is_proper
		end

		@testset "Maps Structure" begin
			# Test Maps structure
			map = Maps(1, 2)
			@test map.atom == 1
			@test map.translation == 2
		end

		@testset "Symmetry Operation Comparison" begin
			# Test isless for SymmetryOperation
			symop1 = SymmetryOperation(
				SMatrix{3, 3, Float64}([
					1.0 0.0 0.0;
					0.0 1.0 0.0;
					0.0 0.0 1.0
				]),  # rotation_frac
				SMatrix{3, 3, Float64}([
					1.0 0.0 0.0;
					0.0 1.0 0.0;
					0.0 0.0 1.0
				]),  # rotation_cart
				SVector{3, Float64}([0.0, 0.0, 0.0]),  # translation_frac
				false, false, true,
			)

			symop2 = SymmetryOperation(
				SMatrix{3, 3, Float64}([
					1.0 0.0 0.0;
					0.0 1.0 0.0;
					0.0 0.0 1.0
				]),  # rotation_frac
				SMatrix{3, 3, Float64}([
					1.0 0.0 0.0;
					0.0 1.0 0.0;
					0.0 0.0 1.0
				]),  # rotation_cart
				SVector{3, Float64}([0.1, 0.0, 0.0]),  # translation_frac
				false, true, true,
			)

			@test symop1 < symop2
		end
	end

	@testset "Simple Structures" begin
		@testset "1x1x1 Simple Cubic (Polonium)" begin
			# Create a simple cubic structure
			lattice_vectors = SMatrix{3, 3, Float64}([
				1.0 0.0 0.0;
				0.0 1.0 0.0;
				0.0 0.0 1.0
			])
			is_periodic = SVector{3, Bool}([true, true, true])
			kd_name = ["Po"]
			kd_int_list = [1]
			x_frac = reshape([0.0, 0.0, 0.0], 3, 1)

			structure = Structure(
				lattice_vectors,
				is_periodic,
				kd_name,
				kd_int_list,
				x_frac,
				verbosity = false
			)

			# Create symmetry with tolerance
			symmetry = Symmetry(structure, 1e-5, verbosity = false)

			# Test basic properties
			@test symmetry.international_symbol == "Pm-3m"
			@test symmetry.spacegroup_number == 221  # Pm-3m space group for simple cubic
			@test symmetry.nat_prim == 1  # One atom in primitive cell
			@test symmetry.nsym == 48  # Number of symmetry operations for cubic
			@test symmetry.ntran == 1
			@test !isempty(symmetry.symdata)
			@test length(symmetry.atoms_in_prim) == symmetry.nat_prim
		end

		@testset "2x2x2 Simple Cubic (Polonium)" begin
			# Create a 2x2x2 supercell of simple cubic structure
			lattice_vectors = SMatrix{3, 3, Float64}([
				2.0 0.0 0.0;
				0.0 2.0 0.0;
				0.0 0.0 2.0
			])
			is_periodic = SVector{3, Bool}([true, true, true])
			kd_name = ["Po" for _ in 1:8]  # 8 atoms in the supercell
			kd_int_list = ones(Int, 8)     # All atoms are Po (type 1)
			
			# Create fractional coordinates for 2x2x2 supercell
			x_frac = [
				0.0 0.0 0.0 0.0 0.5 0.5 0.5 0.5;
				0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.5;
				0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5
			]

			structure = Structure(
				lattice_vectors,
				is_periodic,
				kd_name,
				kd_int_list,
				x_frac,
				verbosity = false
			)

			# Create symmetry with tolerance
			symmetry = Symmetry(structure, 1e-5, verbosity = false)

			# Test basic properties
			@test symmetry.international_symbol == "Pm-3m"
			@test symmetry.spacegroup_number == 221  # Pm-3m space group for simple cubic
			@test symmetry.nat_prim == 1  # Still one atom in primitive cell
			@test symmetry.nsym == 384    # Number of symmetry operations for supercell
			@test symmetry.ntran == 8    # Number of translations for 2x2x2 supercell
			@test !isempty(symmetry.symdata)
			@test length(symmetry.atoms_in_prim) == symmetry.nat_prim

			# Test that atoms are correctly mapped by identity operation
			for i in 1:8
				@test symmetry.symdata[1].rotation_frac * SVector{3, Float64}(x_frac[:, i]) +
					  symmetry.symdata[1].translation_frac ≈ SVector{3, Float64}(x_frac[:, i]) atol=1e-5
			end

		end

		@testset "Two-Atom Structure (CsCl-type)" begin
			# Create a CsCl-type structure
			lattice_vectors = SMatrix{3, 3, Float64}([
				1.0 0.0 0.0;
				0.0 1.0 0.0;
				0.0 0.0 1.0
			])
			is_periodic = SVector{3, Bool}([true, true, true])
			kd_name = ["Cs", "Cl"]
			kd_int_list = [1, 2]
			x_frac = [0.0 0.5;
					  0.0 0.5;
					  0.0 0.5]

			structure = Structure(
				lattice_vectors,
				is_periodic,
				kd_name,
				kd_int_list,
				x_frac,
				verbosity = false
			)

			# Create symmetry with tolerance
			symmetry = Symmetry(structure, 1e-5, verbosity = false)

			# Test basic properties
			@test symmetry.international_symbol == "Pm-3m"
			@test symmetry.spacegroup_number == 221  # Pm-3m space group for CsCl structure
			@test symmetry.nat_prim == 2  # Two atoms in primitive cell
			@test symmetry.nsym == 48  # Number of symmetry operations for cubic
			@test symmetry.ntran == 1
			@test !isempty(symmetry.symdata)
			@test length(symmetry.atoms_in_prim) == symmetry.nat_prim

			# Test that both atoms are correctly mapped by identity operation
			@test symmetry.symdata[1].rotation_frac * SVector{3, Float64}(x_frac[:, 1]) +
				  symmetry.symdata[1].translation_frac ≈ SVector{3, Float64}(x_frac[:, 1]) atol=1e-5

			@test symmetry.symdata[1].rotation_frac * SVector{3, Float64}(x_frac[:, 2]) +
				  symmetry.symdata[1].translation_frac ≈ SVector{3, Float64}(x_frac[:, 2]) atol=1e-5
		end
	end

	@testset "Error Cases" begin
		@testset "Invalid Symmetry Construction" begin
			# Test with invalid tolerance
			lattice_vectors = SMatrix{3, 3, Float64}([
				1.0 0.0 0.0;
				0.0 1.0 0.0;
				0.0 0.0 1.0
			])
			is_periodic = SVector{3, Bool}([true, true, true])
			kd_name = ["Fe"]
			kd_int_list = [1]
			x_frac = reshape([0.0, 0.0, 0.0], 3, 1)

			structure = Structure(
				lattice_vectors,
				is_periodic,
				kd_name,
				kd_int_list,
				x_frac,
				verbosity = false
			)

			@test_throws ArgumentError Symmetry(structure, -1.0, verbosity = false)
		end
	end
end
