include("../src/common/SortedContainer.jl")
include("../src/types/AtomCells.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/utils/RotationMatrices.jl")

using Magesty
using Test

@testset "Component tests" begin
	include("./test_SortedContainer.jl")
	include("./test_AtomicIndices.jl")
	include("./test_UnitaryMatrixCl.jl")
	include("./test_RotationMatrices.jl")
end
# @testset "Magesty.jl" begin
# end
