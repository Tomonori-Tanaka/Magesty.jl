using TOML
using Printf
using Test
using StaticArrays
using Magesty.Structures
using Magesty.Symmetries


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
input["regression"]["datafile"] = joinpath(@__DIR__, "EMBSET.dat")

system = System(input, verbosity = false)

const NUM_CELLS = 27  # Total number of cells: center cell and its neighboring virtual cells
@testset "Fe BCC 2x2x2 Tests" begin
	@testset "Structure Tests" begin
		structure = system.structure

		# Test lattice and basic properties
		a = 2.83  # Fe lattice constant in angstrom
		expected_lattice = SMatrix{3, 3, Float64}([
			2*a  0.0  0.0;
			0.0  2*a  0.0;
			0.0  0.0  2*a
		])
		@test structure.supercell.lattice_vectors ≈ expected_lattice atol=1e-5

		expected_reciprocal = SMatrix{3, 3, Float64}([
			1/(2*a)  0.0  0.0;
			0.0  1/(2*a)  0.0;
			0.0  0.0  1/(2*a)
		])
		@test structure.supercell.reciprocal_vectors ≈ expected_reciprocal atol=1e-5


		# Test basic properties
		@test structure.supercell.num_atoms == 16  # 2x2x2 supercell of BCC (2 atoms/unit cell)
		@test structure.supercell.num_elements == 1  # Only Fe atoms
		@test structure.supercell.kd_int_list == [1 for _ in 1:structure.supercell.num_atoms]
		@test all(structure.is_periodic)  # Periodic in all directions

		# Test atomic positions
		# BCC positions in 2x2x2 supercell
		expected_positions_vvec = [
			[0.00, 0.00, 0.00],
			[0.25, 0.25, 0.25],
			[0.00, 0.00, 0.50],
			[0.50, 0.00, 0.00],
			[0.00, 0.50, 0.00],
			[0.25, 0.25, 0.75],
			[0.75, 0.25, 0.25],
			[0.25, 0.75, 0.25],
			[0.00, 0.50, 0.50],
			[0.50, 0.00, 0.50],
			[0.50, 0.50, 0.00],
			[0.25, 0.75, 0.75],
			[0.75, 0.25, 0.75],
			[0.75, 0.75, 0.25],
			[0.50, 0.50, 0.50],
			[0.75, 0.75, 0.75],
		]
		expected_positions = reduce(hcat, expected_positions_vvec)
		@test structure.supercell.x_frac ≈ expected_positions atol=1e-5

		@test structure.is_periodic == [true, true, true]
		@test structure.kd_name == ["Fe"]
		@test structure.x_image_frac[:, :, 1] ≈ structure.supercell.x_frac atol=1e-5
		@test structure.x_image_frac[:, :, 2] ≈ structure.supercell.x_frac .- 1.0 atol=1e-5
		@test structure.x_image_frac[:, :, NUM_CELLS] ≈ structure.supercell.x_frac .+ 1.0 atol=1e-5
		@test structure.x_image_cart[:, :, 1] ≈
			  structure.supercell.lattice_vectors * structure.supercell.x_frac atol=1e-5
		@test structure.x_image_cart[:, :, 2] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .- 1.0) atol=1e-5
		@test structure.x_image_cart[:, :, NUM_CELLS] ≈
			  structure.supercell.lattice_vectors * (structure.supercell.x_frac .+ 1.0) atol=1e-5
		@test structure.exist_image == [true for _ in 1:NUM_CELLS]
		@test structure.atomtype_group ==
			  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]
	end

	@testset "Symmetry Tests" begin
		symmetry = system.symmetry

		# Test basic properties
		@test symmetry.international_symbol == "Im-3m"  # Space group for BCC Fe
		@test symmetry.spacegroup_number == 229
		@test symmetry.nat_prim == 1  # Two atoms in primitive cell
		@test symmetry.nsym == 768  # 48 point group operations × 16 translations for 2x2x2
		@test symmetry.ntran == 16   # Number of translations for 2x2x2 supercell
		@test !isempty(symmetry.symdata)
		@test length(symmetry.atoms_in_prim) == symmetry.nat_prim

		# Test that identity operation preserves atomic positions
		structure = system.structure
		identity_op = symmetry.symdata[1]
		for i in 1:structure.supercell.num_atoms
			pos = SVector{3, Float64}(structure.supercell.x_frac[:, i])
			transformed_pos = identity_op.rotation_frac * pos + identity_op.translation_frac
			@test transformed_pos ≈ pos atol=1e-5
		end
	end

	sclus = SpinCluster(system, input, verbosity = false)
	Magesty.write_xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))
	structure = Structure(joinpath(@__DIR__, "scecoeffs.xml"), verbosity = false)

	@testset "calc_energy" begin
		spin_config_list = sclus.optimize.spinconfig_list
		energy_list_from_salc::Vector{Float64} = Vector{Float64}(undef, length(spin_config_list))
		for i in eachindex(spin_config_list)
			spin_directions::Matrix{Float64} = spin_config_list[i].spin_directions
			energy_list_from_salc[i] = Magesty.calc_energy(sclus, spin_directions)
			@test abs(energy_list_from_salc[i] - spin_config_list[i].energy) < 0.1
		end
	end
end
