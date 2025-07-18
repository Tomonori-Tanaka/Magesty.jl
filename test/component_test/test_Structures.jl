using Test
using StaticArrays
using Magesty.Structures

@testset "Structure Tests" begin
    @testset "Cell Construction" begin
        # Test for basic cubic lattice
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        kd_int_list = [1, 1]  # Two atoms of the same type
        x_frac = [
            0.0 0.5;
            0.0 0.5;
            0.0 0.0
        ]
        
        cell = Magesty.Structures.Cell(lattice_vectors, kd_int_list, x_frac)
        
        @test cell.num_atoms == 2
        @test cell.num_elements == 1
        @test cell.lattice_vectors == lattice_vectors
        @test cell.reciprocal_vectors ≈ inv(lattice_vectors)
    end

    @testset "Structure Construction" begin
        # Test for basic structure
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        is_periodic = SVector{3,Bool}([true, true, true])
        kd_name = ["Fe"]
        kd_int_list = [1, 1]
        x_frac = [
            0.0 0.5;
            0.0 0.5;
            0.0 0.0
        ]

        structure = Structure(
            lattice_vectors,
            is_periodic,
            kd_name,
            kd_int_list,
            x_frac,
            verbosity = false
        )

        @test structure.supercell.num_atoms == 2
        @test structure.supercell.num_elements == 1
        @test structure.is_periodic == is_periodic
        @test structure.kd_name == kd_name
        @test size(structure.x_image_frac) == (3, 2, 27)  # 3 dimensions × 2 atoms × 27 cells
        @test size(structure.x_image_cart) == (3, 2, 27)
    end

    @testset "Invalid Lattice Vectors" begin
        # Test for linearly dependent lattice vectors
        invalid_lattice = SMatrix{3,3,Float64}([
            1.0 0.0 1.0;
            0.0 1.0 0.0;
            0.0 0.0 0.0
        ])
        kd_int_list = [1]
        x_frac = [0.0; 0.0; 0.0;;]

        @test_throws ErrorException Structures.Cell(invalid_lattice, kd_int_list, x_frac)
    end

    @testset "Atom Type Grouping" begin
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        is_periodic = SVector{3,Bool}([true, true, true])
        kd_name = ["Fe", "Co"]
        kd_int_list = [1, 2, 1, 2]  # Fe-Co-Fe-Co pattern
        x_frac = [
            0.0 0.25 0.5 0.75;
            0.0 0.0  0.0 0.0;
            0.0 0.0  0.0 0.0
        ]

        structure = Structure(
            lattice_vectors,
            is_periodic,
            kd_name,
            kd_int_list,
            x_frac,
            verbosity = false
        )

        # Test atom type grouping
        @test length(structure.atomtype_group) == 2  # Two types: Fe and Co
        @test structure.atomtype_group[1] == [1, 3]  # Fe positions
        @test structure.atomtype_group[2] == [2, 4]  # Co positions
    end
end 