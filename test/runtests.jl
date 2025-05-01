include("../src/common/SortedContainer.jl")
include("../src/common/CountingContainer.jl")
include("../src/types/AtomCell.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/types/SpinConfig.jl")
include("../src/utils/RotationMatrix.jl")
include("../src/utils/MySphericalHarmonics.jl")

using Magesty
using Test

@testset "component tests" begin
	include("./component_test/test_MySphericalHarmonics.jl")
	include("./component_test/test_SortedContainer.jl")
	include("./component_test/test_CountingContainer.jl")
	include("./component_test/test_AtomicIndices.jl")
	include("./component_test/test_UnitaryMatrixCl.jl")
	include("./component_test/test_RotationMatrix.jl")
	include("./component_test/test_SpinConfig.jl")
end

@testset "examples" begin
	include("./examples/febcc_2x2x2_pm/test.jl")
	# include("./examples/fept_tetragonal_2x2x2/test.jl")
	# include("./examples/fecob2_3x3x3/test.jl")
end
