include("../src/common/SortedContainer.jl")
include("../src/common/CountingContainer.jl")
include("../src/types/AtomCell.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/utils/RotationMatrix.jl")

using Magesty
using Test

@testset "component tests" begin
	include("./component_test/test_SortedContainer.jl")
	include("./component_test/test_CountingContainer.jl")
	include("./component_test/test_AtomicIndices.jl")
	include("./component_test/test_UnitaryMatrixCl.jl")
	include("./component_test/test_RotationMatrix.jl")
end

@testset "examples" begin
	# include("./examples/b2feco1x1x1.jl")
	# include("./examples/bccfe1x1x1.jl")
	# include("./examples/b2feco1x1x2.jl")
	include("./examples/bccfe2x2x2.jl")
	# include("./examples/b2feco2x2x2.jl")
	# include("./examples/fesi_atomicchain.jl")
	# include("./examples/hcpco1x1x1.jl")
end
# @testset "Magesty.jl" begin
# end
