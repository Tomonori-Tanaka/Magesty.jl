using Test
using Magesty
using TOML

@testset "FePt L1_0 2x2x2 Tests" begin
	input = TOML.parsefile(joinpath(@__DIR__, "input.toml"))

	# Build, save, and reload an SCEBasis. The reload path must reproduce the
	# SALC structure exactly so downstream design matrices match.
	basis = SCEBasis(input; verbosity = false)

	system_xml = joinpath(@__DIR__, "system.xml")
	save(basis, system_xml)
	basis_from_xml = load(SCEBasis, system_xml)

	@testset "SCEBasis XML round-trip" begin
		@test length(basis.salcbasis.salc_list) ==
			  length(basis_from_xml.salcbasis.salc_list)
		for (i, (orig, xml)) in enumerate(
			zip(basis.salcbasis.salc_list, basis_from_xml.salcbasis.salc_list),
		)
			@test length(orig) == length(xml)
			for (j, (cbc_orig, cbc_xml)) in enumerate(zip(orig, xml))
				@test cbc_orig.ls == cbc_xml.ls
				@test cbc_orig.Lf == cbc_xml.Lf
				@test cbc_orig.Lseq == cbc_xml.Lseq
				@test cbc_orig.atoms == cbc_xml.atoms
				@test cbc_orig.multiplicity == cbc_xml.multiplicity
				@test isapprox(cbc_orig.coefficient, cbc_xml.coefficient, atol = 1e-10)
				@test size(cbc_orig.coeff_tensor) == size(cbc_xml.coeff_tensor)
				@test isapprox(cbc_orig.coeff_tensor, cbc_xml.coeff_tensor, atol = 1e-10)
			end
		end
	end

	# Fit SCE coefficients and persist the fitted model alongside the basis.
	embset_path = joinpath(@__DIR__, "EMBSET.dat")
	dataset = SCEDataset(basis, embset_path)
	fitted = fit(SCEFit, dataset, Ridge(lambda = 0.0); torque_weight = 1.0)
	model = SCEModel(fitted)
	save(model, joinpath(@__DIR__, "scecoeffs.xml"))

	@testset "SCEFit smoke checks" begin
		@test isfinite(intercept(fitted))
		@test length(coef(fitted)) == length(basis.salcbasis.salc_list)
		@test all(isfinite, coef(fitted))
	end
end
