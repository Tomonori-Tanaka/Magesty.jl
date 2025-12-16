using Magesty
using TOML

input = TOML.parse(
	open("./examples/fept_tetragonal_2x2x2/input.toml", "r"),
)

# Build SCE basis using build_sce_basis
system = build_sce_basis(input, verbosity = false)

# Read spin configurations
datafile = joinpath(@__DIR__, "EMBSET.dat")
spinconfig_list = read_embset(datafile, system.structure.supercell.num_atoms)

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

Magesty.write_xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))
