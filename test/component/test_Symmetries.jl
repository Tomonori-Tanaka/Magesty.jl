using Test
using Logging
using StaticArrays
using Magesty.Symmetries
using Magesty.Structures
using Magesty.InputSpecs: SymmetryOptions

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
                true,   # is_proper
            )

            @test symop.rotation_frac == rotation_frac
            @test symop.rotation_cart == rotation_cart
            @test symop.translation_frac == translation_frac
            @test !symop.is_translation
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
                false, true,
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
                false, true,
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

    @testset "Tolerance Drives the spglib Search" begin
        # A cubic CsCl-type cell whose Cl site is displaced along z by a
        # distance that straddles the two tolerances under test.
        #
        #   Cs at (0, 0, 0), Cl at (1/2, 1/2, 1/2 + delta), a = 4.0 Ang.
        #
        # With delta = 1.25e-5 in fractional units the Cl atom sits
        # d = 1.25e-5 * 4.0 = 5.0e-5 Ang off the body center, so any operation
        # that reverses z (the horizontal mirror, the inversion, the two-fold
        # axes normal to z) misplaces it by 2d = 1.0e-4 Ang. spglib accepts an
        # operation when that mismatch is below `symprec`, hence
        #
        #   symprec = 1e-5  ->  1.0e-4 > symprec: z-reversing operations are
        #                       rejected. What survives is the point group of a
        #                       tetragonal polar axis along z, 4mm (order 8:
        #                       E, 2 C4z, C2z, 2 sigma_v, 2 sigma_d), i.e. the
        #                       symmorphic space group P4mm (99).
        #   symprec = 1e-3  ->  1.0e-4 < symprec: the displacement is absorbed
        #                       and the full cubic m-3m (order 48) is recovered,
        #                       i.e. Pm-3m (221).
        #
        # Neither cell is primitive-reducible and both leave two atoms in the
        # primitive cell, so ntran = 1 and nsym equals the point-group order.
        lattice_vectors = SMatrix{3, 3, Float64}([
            4.0 0.0 0.0;
            0.0 4.0 0.0;
            0.0 0.0 4.0
        ])
        is_periodic = SVector{3, Bool}([true, true, true])
        kd_name = ["Cs", "Cl"]
        kd_int_list = [1, 2]
        delta = 1.25e-5
        x_frac = [0.0 0.5;
                  0.0 0.5;
                  0.0 0.5+delta]

        structure = Structure(
            lattice_vectors,
            is_periodic,
            kd_name,
            kd_int_list,
            x_frac,
            verbosity = false
        )

        symmetry_tight = Symmetry(structure, 1e-5, verbosity = false)
        @test symmetry_tight.international_symbol == "P4mm"
        @test symmetry_tight.spacegroup_number == 99
        @test symmetry_tight.nsym == 8
        @test symmetry_tight.ntran == 1
        @test symmetry_tight.nat_prim == 2
        @test symmetry_tight.tol == 1e-5

        symmetry_loose = Symmetry(structure, 1e-3, verbosity = false)
        @test symmetry_loose.international_symbol == "Pm-3m"
        @test symmetry_loose.spacegroup_number == 221
        @test symmetry_loose.nsym == 48
        @test symmetry_loose.ntran == 1
        @test symmetry_loose.nat_prim == 2
        @test symmetry_loose.tol == 1e-3

        # The `SymmetryOptions` constructor must forward the same tolerance.
        options = SymmetryOptions(tolerance_sym = 1e-3)
        symmetry_from_options = Symmetry(structure, options, verbosity = false)
        @test symmetry_from_options.spacegroup_number ==
              symmetry_loose.spacegroup_number
        @test symmetry_from_options.nsym == symmetry_loose.nsym
    end

    @testset "Looser-Tolerance Warning" begin
        # Same cell as above: a = 4.0 Ang, Cl displaced along z by
        # d = 1.25e-5 * 4.0 = 5.0e-5 Ang, so a z-reversing operation misplaces
        # it by 2d = 1.0e-4 Ang. The re-search runs at ten times the requested
        # tolerance, which brackets that mismatch for tol = 3e-5:
        #
        #   tol      = 3.0e-5 < 1.0e-4  -> P4mm  (99),  8 operations
        #   10 * tol = 3.0e-4 > 1.0e-4  -> Pm-3m (221), 48 operations
        #
        # so the looser search strictly gains operations and must warn. At
        # tol = 1e-3 both searches already sit above the mismatch, return the
        # same 48 operations, and must stay silent.
        lattice_vectors = SMatrix{3, 3, Float64}([
            4.0 0.0 0.0;
            0.0 4.0 0.0;
            0.0 0.0 4.0
        ])
        is_periodic = SVector{3, Bool}([true, true, true])
        kd_name = ["Cs", "Cl"]
        kd_int_list = [1, 2]
        delta = 1.25e-5
        x_frac = [0.0 0.5;
                  0.0 0.5;
                  0.0 0.5+delta]

        structure = Structure(
            lattice_vectors,
            is_periodic,
            kd_name,
            kd_int_list,
            x_frac,
            verbosity = false
        )

        # The warning is a `verbosity` diagnostic; the extra search is skipped
        # entirely when the caller asked for silence.
        redirect_stdout(devnull) do
            @test_logs (:warn,) match_mode = :any Symmetry(
                structure, 3e-5, verbosity = true)
            @test_logs min_level = Logging.Warn Symmetry(
                structure, 1e-3, verbosity = true)
        end
        @test_logs min_level = Logging.Warn Symmetry(
            structure, 3e-5, verbosity = false)
    end

    @testset "Symmetry Printout Reports the Tolerance" begin
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

        printout = mktemp() do path, io
            redirect_stdout(io) do
                Symmetry(structure, 2.5e-4, verbosity = true)
            end
            flush(io)
            return read(path, String)
        end
        @test occursin("Symmetry tolerance (symprec) = 2.500e-04", printout)
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
