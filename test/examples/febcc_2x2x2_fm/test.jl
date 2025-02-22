using TOML

input = TOML.parse(
	open("/Users/tomorin/Packages/Magesty/test/examples/febcc_2x2x2_fm/input.toml", "r"),
)
sclus = SpinCluster(input)
@show sclus.cluster.cluster_list_with_cell
@show sclus.basisset.basislist
display(sclus.basisset.classified_basisdict)
display(sclus.basisset.salc_list)
display(sclus.optimize.SCE)