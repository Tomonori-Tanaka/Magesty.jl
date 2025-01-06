include("../src/common/SortedContainer.jl")
include("../src/types/AtomCell.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/utils/RotationMatrix.jl")

using Magesty
using Test

@testset "component tests" begin
	include("./component_test/test_SortedContainer.jl")
	include("./component_test/test_AtomicIndices.jl")
	include("./component_test/test_UnitaryMatrixCl.jl")
	include("./component_test/test_RotationMatrix.jl")
end

@testset "examples" begin
	include("./examples/b2feco1x1x1.jl")
end
# @testset "Magesty.jl" begin
# end
