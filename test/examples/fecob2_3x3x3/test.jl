using TOML

input = TOML.parse(open("/Users/tomorin/Packages/Magesty/test/examples/fecob2_3x3x3/input.toml", "r"))
sclus = SpinCluster(input)

@show sclus.cluster.cluster_list
@show sclus.basisset.basislist
display(sclus.basisset.classified_basisdict)
display(sclus.basisset.salc_list)
display(sclus.optimize.SCE)