include("../src/common/SortedContainer.jl")
include("../src/types/AtomicIndices.jl")

using Magesty
using Test

include("./test_SortedContainer.jl")
include("./test_AtomicIndices.jl")

@testset "Magesty.jl" begin
    # Write your tests here.
end
