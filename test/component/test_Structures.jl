using Test
using LinearAlgebra: det, I
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

    # ------------------------------------------------------------------
    # Coverage uplift: surface area below was previously untested.
    # ------------------------------------------------------------------

    @testset "Left-handed lattice rejected" begin
        # _validate_lattice_vectors rejects det < 0 (not just det ≈ 0).
        # Swapping two columns of the identity flips the orientation.
        left_handed = SMatrix{3,3,Float64}([
            0.0 1.0 0.0;
            1.0 0.0 0.0;
            0.0 0.0 1.0
        ])  # det = -1
        @test det(left_handed) < 0
        kd_int_list = [1]
        x_frac = [0.0; 0.0; 0.0;;]
        @test_throws ErrorException Magesty.Structures.Cell(
            left_handed, kd_int_list, x_frac,
        )
    end

    @testset "Cell show" begin
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        cell = Magesty.Structures.Cell(lattice_vectors, [1, 1], [0.0 0.5; 0.0 0.5; 0.0 0.0])
        s = sprint(show, cell)
        # Each field documented in the source's show body appears.
        @test occursin("lattice_vectors:", s)
        @test occursin("reciprocal_vectors:", s)
        @test occursin("num_atoms:", s)
        @test occursin("num_elements:", s)
        @test occursin("kd_ind_list:", s)
        @test occursin("x_frac:", s)
    end

    @testset "Structure show" begin
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        structure = Structure(
            lattice_vectors,
            SVector{3,Bool}([true, true, true]),
            ["Fe"],
            [1, 1],
            [0.0 0.5; 0.0 0.5; 0.0 0.0],
            verbosity = false,
        )
        s = sprint(show, structure)
        @test occursin("supercell:", s)
        @test occursin("is_periodic:", s)
        @test occursin("kd_name:", s)
        @test occursin("x_image_frac:", s)
        @test occursin("x_image_cart:", s)
        @test occursin("exist_image:", s)
        @test occursin("atomtype_group:", s)
    end

    @testset "verbosity=true emits structure info" begin
        # The verbose path calls _print_structure_stdout + a "Time Elapsed"
        # line. Capture via a Pipe — IOBuffer is not supported by
        # redirect_stdout on Julia 1.12 (needs a file descriptor).
        lattice_vectors = SMatrix{3,3,Float64}([
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])
        orig_stdout = stdout
        rd, wr = redirect_stdout()
        try
            Structure(
                lattice_vectors,
                SVector{3,Bool}([true, true, true]),
                ["Fe", "Co"],
                [1, 2],
                [0.0 0.5; 0.0 0.5; 0.0 0.0],
                verbosity = true,
            )
        finally
            close(wr)
            redirect_stdout(orig_stdout)
        end
        captured = read(rd, String)
        @test occursin("SYSTEM", captured)
        @test occursin("Lattice vector", captured)
        @test occursin("Atomic species:", captured)
        @test occursin("Atomic positions in fractional basis:", captured)
        @test occursin("Time Elapsed:", captured)
        @test occursin("Fe", captured)
        @test occursin("Co", captured)
    end

    @testset "Structure from XML" begin
        # Minimal valid XML matching the schema the constructor walks.
        valid_xml = """<?xml version="1.0"?>
        <root>
          <System>
            <LatticeVector>
              <a1>1.0 0.0 0.0</a1>
              <a2>0.0 1.0 0.0</a2>
              <a3>0.0 0.0 1.0</a3>
            </LatticeVector>
            <Periodicity>1 1 1</Periodicity>
            <Positions>
              <pos element="Fe">0.0 0.0 0.0</pos>
              <pos element="Fe">0.5 0.5 0.5</pos>
            </Positions>
          </System>
        </root>
        """

        function write_to_temp(content::AbstractString)
            path, io = mktemp()
            write(io, content)
            close(io)
            return path
        end

        @testset "Happy path" begin
            path = write_to_temp(valid_xml)
            structure = Structure(path; verbosity = false)
            @test structure.supercell.num_atoms == 2
            @test structure.supercell.num_elements == 1
            @test structure.kd_name == ["Fe"]
            @test structure.is_periodic == SVector{3, Bool}([true, true, true])
            @test structure.supercell.lattice_vectors == SMatrix{3, 3, Float64}(I)
        end

        @testset "Mixed elements sorted alphabetically" begin
            # kd_name must be the unique element set sorted (per source
            # L274). Co < Fe alphabetically, so kd_name == ["Co", "Fe"].
            xml = replace(valid_xml,
                "<pos element=\"Fe\">0.5 0.5 0.5</pos>" =>
                "<pos element=\"Co\">0.5 0.5 0.5</pos>")
            path = write_to_temp(xml)
            structure = Structure(path; verbosity = false)
            @test structure.kd_name == ["Co", "Fe"]
            @test structure.supercell.num_elements == 2
            # First atom is Fe → index 2 in sorted kd_name. Second is Co
            # → index 1.
            @test structure.supercell.kd_int_list == [2, 1]
        end

        # --- error paths: one per throw site in the constructor body ---

        @testset "Missing <System> node" begin
            xml = replace(valid_xml, r"<System>.*</System>"s => "")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Missing <LatticeVector> node" begin
            xml = replace(valid_xml, r"<LatticeVector>.*</LatticeVector>"s => "")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Missing <a1> node" begin
            xml = replace(valid_xml, "<a1>1.0 0.0 0.0</a1>" => "")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Invalid lattice vector format" begin
            # Wrong number of components in a vector.
            xml = replace(valid_xml,
                "<a1>1.0 0.0 0.0</a1>" => "<a1>1.0 0.0</a1>")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Missing <Periodicity> node" begin
            xml = replace(valid_xml, "<Periodicity>1 1 1</Periodicity>" => "")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Invalid periodicity format" begin
            xml = replace(valid_xml,
                "<Periodicity>1 1 1</Periodicity>" =>
                "<Periodicity>1 1</Periodicity>")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Missing <Positions> node" begin
            xml = replace(valid_xml, r"<Positions>.*</Positions>"s => "")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Empty <Positions>" begin
            xml = replace(valid_xml,
                r"<Positions>.*</Positions>"s =>
                "<Positions></Positions>")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end

        @testset "Invalid position format" begin
            xml = replace(valid_xml,
                "<pos element=\"Fe\">0.5 0.5 0.5</pos>" =>
                "<pos element=\"Fe\">0.5 0.5</pos>")
            path = write_to_temp(xml)
            @test_throws ArgumentError Structure(path; verbosity = false)
        end
    end
end