include("../src/SortedCounters.jl")
include("../src/AtomCells.jl")
include("../src/SpinConfigs.jl")
include("../src/SphericalHarmonicsTransforms.jl")
include("../src/RotationMatrix.jl")
include("../src/TesseralHarmonics.jl")
include("../src/InputSpecs.jl")
include("../src/AngularMomentumCoupling.jl")
include("../src/CoupledBases.jl")
include("helpers/fileutils.jl")

using Magesty
using Test

const TEST_MODE = get(ENV, "TEST_MODE", "all")

if TEST_MODE in ("unit", "all")
	@testset verbose = true "component tests" begin
		include("./component_test/test_TesseralHarmonics.jl")
		include("./component_test/test_AngularMomentumCoupling.jl")
		include("./component_test/test_SortedCounter.jl")
		include("./component_test/test_SphericalHarmonicsTransforms.jl")
		include("./component_test/test_RotationMatrix.jl")
		include("./component_test/test_SpinConfigs.jl")
		include("./component_test/test_InputSpecs.jl")
		include("./component_test/test_cluster_cutoff_criterion.jl")
		include("./component_test/test_Structures.jl")
		include("./component_test/test_Symmetries.jl")
		include("./component_test/test_CoupledBases.jl")
		include("./component_test/test_Fitting_dispatch.jl")
		include("./component_test/test_SALCBases_l13_regression.jl")
		include("./component_test/test_salc_canonical_gauge.jl")
		include("./component_test/test_SCEBasis.jl")
		include("./component_test/test_SCEDataset.jl")
		include("./component_test/test_SCEFit.jl")
		include("./component_test/test_save_load.jl")
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
if TEST_MODE in ("sphericart",)
	@testset verbose = true "SpheriCart agreement" begin
		include("./component_test/test_sphericart_agreement.jl")
	end
end
