using Aqua
using Test
using Magesty

@testset "Aqua package quality" begin
    Aqua.test_all(
        Magesty;
        ambiguities = false,  # external packages can introduce ambiguities
    )
end
