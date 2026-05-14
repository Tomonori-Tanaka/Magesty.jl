include("../src/common/SortedCounter.jl")
include("../src/types/AtomCells.jl")
include("../src/SpinConfigs.jl")
include("../src/utils/SphericalHarmonicsTransforms.jl")
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
	@testset verbose = true "component tests" begin
		include("./component_test/test_MySphericalHarmonics.jl")
		include("./component_test/test_AngularMomentumCoupling.jl")
		include("./component_test/test_SortedCounter.jl")
		include("./component_test/test_SphericalHarmonicsTransforms.jl")
		include("./component_test/test_RotationMatrix.jl")
		include("./component_test/test_SpinConfigs.jl")
		include("./component_test/test_ConfigParser.jl")
		include("./component_test/test_Structures.jl")
		include("./component_test/test_Symmetries.jl")
		include("./component_test/test_Basis.jl")
		# include("./component_test/test_Optimize.jl")
		include("./component_test/test_Optimize_dispatch.jl")
		include("./component_test/test_SALCBases_l13_regression.jl")
		include("./component_test/test_SCEBasis.jl")
		include("./component_test/test_Version.jl")
	end
end
if TEST_MODE in ("integration", "all")
	@testset verbose = true "examples" begin
		include("./examples/square_lattice/test.jl")
		include("./examples/dimer/test.jl")
		include("./examples/chain/test.jl")
		include("./examples/febcc_2x2x2_pm/test.jl")
		include("./examples/fept_tetragonal_2x2x2/test.jl")
		include("./examples/fege_2x2x2/test.jl")
	end
end
if TEST_MODE in ("develop",)
	@testset verbose = true "develop" begin
		include("./develop_tmp/test_develop.jl")
	end
end
if TEST_MODE in ("jet",)
	include("./jet.jl")
end
if TEST_MODE in ("aqua",)
	include("./aqua.jl")
end
if TEST_MODE in ("bench_sphericart",)
	include("./benchmark_sphericart.jl")
end
if TEST_MODE in ("bench_optimize_sphericart",)
	include("./benchmark_optimize_sphericart.jl")
end
if TEST_MODE in ("sphericart",)
	@testset verbose = true "SpheriCart agreement" begin
		include("./component_test/test_sphericart_agreement.jl")
	end
end
