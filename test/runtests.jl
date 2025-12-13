include("../src/common/SortedContainer.jl")
include("../src/common/CountingContainer.jl")
include("../src/types/AtomCells.jl")
include("../src/types/AtomicIndices.jl")
include("../src/types/UnitaryMatrixCl.jl")
include("../src/SpinConfigs.jl")
include("../src/utils/RotationMatrix.jl")
include("../src/utils/MySphericalHarmonics.jl")
include("../src/utils/ConfigParser.jl")
include("../src/utils/AngularMomentumCoupling.jl")
include("../src/types/Basis.jl")
include("helpers/fileutils.jl")

using Magesty
using Test

const TEST_MODE = get(ENV, "TEST_MODE", "all")

if TEST_MODE in ("unit", "all")
	@testset "component tests" begin
		include("./component_test/test_MySphericalHarmonics.jl")
		include("./component_test/test_AngularMomentumCoupling.jl")
		include("./component_test/test_SortedContainer.jl")
		include("./component_test/test_CountingContainer.jl")
		include("./component_test/test_AtomicIndices.jl")
		include("./component_test/test_UnitaryMatrixCl.jl")
		include("./component_test/test_RotationMatrix.jl")
		include("./component_test/test_SpinConfigs.jl")
		include("./component_test/test_ConfigParser.jl")
		include("./component_test/test_Structures.jl")
		include("./component_test/test_Symmetries.jl")
		include("./component_test/Basis.jl")
	end
end
if TEST_MODE in ("integration", "all")
	@testset "examples" begin
		include("./examples/dimer/test.jl")
		# include("./examples/febcc_2x2x2_pm/test.jl")
		# include("./examples/fept_tetragonal_2x2x2/test.jl")
	end
end
if TEST_MODE in ("develop",)
	@testset "develop" begin
		include("./develop_tmp/test_develop.jl")
	end
end
