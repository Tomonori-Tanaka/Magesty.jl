using TOML

io = open("./result/fept_tetragonal_2x2x2.txt", "w")

input = TOML.parse(
	open("./examples/fept_tetragonal_2x2x2/input.toml", "r"),
)
system = System(input, verbosity = false)

# print info to file
redirect_stdout(io) do
	print_info(system)
end

sclus = SpinCluster(system, input, verbosity = false)
Magesty.write_sce2xml(sclus, joinpath(@__DIR__, "scecoeffs.xml"))

close(io)