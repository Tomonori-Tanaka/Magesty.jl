include("../src/common/SortedContainer.jl")
include("../src/types/AtomCells.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/utils/RotationMatrices.jl")

using Magesty
using Test

@testset "component tests" begin
	include("./component_test/test_SortedContainer.jl")
	include("./component_test/test_AtomicIndices.jl")
	include("./component_test/test_UnitaryMatrixCl.jl")
	include("./component_test/test_RotationMatrices.jl")
end
# @testset "Magesty.jl" begin
# end
