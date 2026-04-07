using TOML
using Printf
using Test
using LinearAlgebra
using Magesty


@testset "Square Lattice Tests" begin

	spin_directions_fm = [
		0.0 0.0 0.0 0.0;
		0.0 0.0 0.0 0.0;
		1.0 1.0 1.0 1.0]
	spin_directions_afm = [
		0.0  0.0  0.0 0.0;
		0.0  0.0  0.0 0.0;
		1.0 -1.0 -1.0 1.0]

	@testset "Hypothetical SCE model (isotropic)" begin
		input = TOML.parse(open(joinpath(@__DIR__, "input_isotropic.toml"), "r"))
		system = System(input, verbosity = false)
		Magesty.write_xml(system, joinpath(@__DIR__, "system_isotropic.xml"))
		# Heisenberg model: E = Σ_{<ij>} Jij * Si·Sj, Jij = -1.0 eV (ferromagnetic) for all pairs.
		#
		# Square lattice (4 atoms, PBC): 8 NN pairs + 8 diagonal (NNN) pairs within the supercell.
		# Each SALC has a symmetry coefficient c = norm of the first basis's coefficient vector.
		# For Lf=0, the coefficient vector has length 1, so c = coefficient[1].
		# For Lf=2, the vector has length 5 (one per Mf), so we take the norm.
		#
		# Conversion from physical Jij to SCE coefficient jphi:
		#   jphi = Jij / (c × √3)
		# where the √3 comes from the l=1 spherical harmonics coupling:
		#   (4π) × (1/√3) × (3/4π) = √3  per pair.
		J = -1.0
		# Read SALC coefficients from system to be robust against ordering changes.
		# When a SALC has multiple bases, the norm of the first basis's coefficient vector is used.
		salc_coeffs = [norm(salc[1].coefficient) for salc in system.basisset.salc_list]
		jphi_list = [J / (c * sqrt(3)) for c in salc_coeffs]
		model = SCEModel(
			0.0, # reference_energy
			jphi_list,
			system.basisset,
			system.symmetry,
			system.structure.supercell.num_atoms,
		)

		# Ferromagnetic configuration: all spins along +z, Si·Sj = 1 for all pairs.
		# E = Jij × (8 NN + 8 diagonal) × 1 = -1.0 × 16 = -16 eV
		energy_analytic_fm = J * 8 + J * 8  # -16.0 eV
		@test predict_energy(model, spin_directions_fm) ≈ energy_analytic_fm atol = 1e-6

		# Antiferromagnetic configuration: S1=S4=+z, S2=S3=-z.
		# NN pairs:      all Si·Sj = -1  → Σ_{NN} = -8
		# Diagonal pairs: (1,4) and (2,3) are same spin → Si·Sj = +1 → Σ_{diag} = +8
		# E = J × (-8 + 8) = 0 eV
		energy_analytic_afm = 0.0
		@test predict_energy(model, spin_directions_afm) ≈ energy_analytic_afm atol = 1e-6
	end

	@testset "Hypothetical SCE model (anisotropic)" begin

		# Assuing the exchange tensor considering all pairs:
		exchange_tensor_expected = [
			-1/3 0.0 0.0;
			0.0 -1/3 0.0;
			0.0 0.0 2/3]
		isotropic_jij_expected = 0 # tr(exchange_tensor_expected) / 3
		dm_vector_expected = [0.0, 0.0, 0.0]
		gamma_matrix_expected = [
			-1/3 0.0 0.0;
			0.0 -1/3 0.0;
			0.0 0.0 2/3]

		fm_z_aligned = [
			0.0 0.0 0.0 0.0;
			0.0 0.0 0.0 0.0;
			1.0 1.0 1.0 1.0]
		fm_x_aligned = [
			1.0 1.0 1.0 1.0;
			0.0 0.0 0.0 0.0;
			0.0 0.0 0.0 0.0]

		input = TOML.parse(open(joinpath(@__DIR__, "input_anisotropic.toml"), "r"))
		system = System(input, verbosity = false)
		Magesty.write_xml(system, joinpath(@__DIR__, "system_anisotropic.xml"))

		# Only Lf=2 SALCs contribute to the symmetric anisotropic exchange; set others to zero.
		# Use norm of the coefficient vector (length = 2*Lf+1) to get the scalar SALC weight.
		salc_coeffs = [salc[1].Lf == 2 ? norm(salc[1].coefficient) : 0.0 for salc in system.basisset.salc_list]

		model = SCEModel(
			0.0, # reference_energy
			jphi_list,
			system.basisset,
			system.symmetry,
			system.structure.supercell.num_atoms,
		)


	end

end
