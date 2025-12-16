using Magesty
using TOML

input = TOML.parse(
	open("./examples/fept_tetragonal_2x2x2/input.toml", "r"),
)

# Build SCE basis using build_sce_basis
system = build_sce_basis(input, verbosity = false)

# Write System to XML file
system_xml_file = joinpath(@__DIR__, "system.xml")
write_xml(system, system_xml_file)

# Test: Compare basisset from XML with original basisset
system_from_xml = build_sce_basis_from_xml(input, system_xml_file, verbosity = false)

# Compare salc_list
@assert length(system.basisset.salc_list) == length(system_from_xml.basisset.salc_list) "SALC list length mismatch"

for (i, (salc_orig, salc_xml)) in enumerate(zip(system.basisset.salc_list, system_from_xml.basisset.salc_list))
	@assert length(salc_orig) == length(salc_xml) "SALC group $i length mismatch"
	
	for (j, (cbc_orig, cbc_xml)) in enumerate(zip(salc_orig, salc_xml))
		# Compare basic fields
		@assert cbc_orig.ls == cbc_xml.ls "SALC[$i][$j]: ls mismatch"
		@assert cbc_orig.Lf == cbc_xml.Lf "SALC[$i][$j]: Lf mismatch"
		@assert cbc_orig.Lseq == cbc_xml.Lseq "SALC[$i][$j]: Lseq mismatch"
		@assert cbc_orig.atoms == cbc_xml.atoms "SALC[$i][$j]: atoms mismatch"
		@assert cbc_orig.multiplicity == cbc_xml.multiplicity "SALC[$i][$j]: multiplicity mismatch"
		
		# Compare coefficient vector (numerical)
		@assert isapprox(cbc_orig.coefficient, cbc_xml.coefficient, atol=1e-10) "SALC[$i][$j]: coefficient mismatch"
		
		# Compare coeff_tensor (numerical)
		@assert size(cbc_orig.coeff_tensor) == size(cbc_xml.coeff_tensor) "SALC[$i][$j]: coeff_tensor size mismatch"
		@assert isapprox(cbc_orig.coeff_tensor, cbc_xml.coeff_tensor, atol=1e-10) "SALC[$i][$j]: coeff_tensor mismatch"
	end
end

# Read spin configurations
datafile = joinpath(@__DIR__, "EMBSET.dat")
spinconfig_list = read_embset(datafile)

# Fit SCE coefficients using fit_sce_model
regression_config = input["regression"]
estimator = ElasticNet(alpha=regression_config["alpha"], lambda=regression_config["lambda"])
optimizer = fit_sce_model(
	system,
	spinconfig_list,
	estimator,
	regression_config["weight"],
	verbosity = false,
)

# Create SpinCluster from System and Optimizer
sclus = SpinCluster(
	system.structure,
	system.symmetry,
	system.cluster,
	system.basisset,
	optimizer,
)

xml_file = joinpath(@__DIR__, "scecoeffs.xml")
Magesty.write_xml(sclus, xml_file)

