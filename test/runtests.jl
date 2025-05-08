include("../src/common/SortedContainer.jl")
include("../src/common/CountingContainer.jl")
include("../src/types/AtomCell.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/SpinConfig.jl")
include("../src/utils/RotationMatrix.jl")
include("../src/utils/MySphericalHarmonics.jl")
include("../src/utils/ConfigParser.jl")

using Magesty
using Test

test_mode = get(ENV, "TEST_MODE", "all")

if test_mode in ("unit", "all")
	@testset "component tests" begin
		include("./component_test/test_MySphericalHarmonics.jl")
		include("./component_test/test_SortedContainer.jl")
		include("./component_test/test_CountingContainer.jl")
		include("./component_test/test_AtomicIndices.jl")
		include("./component_test/test_UnitaryMatrixCl.jl")
		include("./component_test/test_RotationMatrix.jl")
		include("./component_test/test_SpinConfig.jl")
		include("./component_test/test_ConfigParser.jl")
		include("./component_test/test_Structure.jl")
	end

	if test_mode in ("integration", "all")
		@testset "examples" begin
			include("./examples/febcc_2x2x2_pm/test.jl")

		end
	end
end

if test_mode in ("develop", "all")
	@testset "develop" begin
		include("./develop_tmp/test_develop.jl")
	end
end
