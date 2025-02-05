using .SpinConfigs
using LinearAlgebra

@testset "SpinConfig" begin
	spinconfigs::Vector{SpinConfig} = read_embset("./examples/feptL10_2x2x2_saxis001/EMBSET.txt", 16)

	@test length(spinconfigs) == 40
	@test isapprox(spinconfigs[1].energy, -121.20626)
	@test isapprox(spinconfigs[2].energy, -121.29061)
	@test isapprox(spinconfigs[40].energy, -121.22161)

	@test isapprox(spinconfigs[1].magmom_size[1], norm([-0.5870000, -0.1960000, 2.0170000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 1], [-0.5870000, -0.1960000, 2.0170000] / norm([-0.5870000, -0.1960000, 2.0170000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 1], [6.20150e-02, 2.55300e-02, 2.07920e-02], atol=1e-4)

    @test isapprox(spinconfigs[1].magmom_size[2], norm([0.6180000, 0.1020000, 2.0140000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 2], [0.6180000, 0.1020000, 2.0140000] / norm([0.6180000, 0.1020000, 2.0140000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 2], [-3.86260e-02, -2.55540e-02, 1.32630e-02], atol=1e-4)

    @test isapprox(spinconfigs[1].magmom_size[16], norm([0.0440000, 0.0120000, 0.2320000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 16], [0.0440000, 0.0120000, 0.2320000] / norm([0.0440000, 0.0120000, 0.2320000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 16], [-1.61470e-01, -4.43760e-02, 3.43780e-02], atol=1e-4)

    @test isapprox(spinconfigs[2].magmom_size[1], norm([0.2640000, 0.4870000, 2.0390000]), atol=1e-4)
    @test isapprox(spinconfigs[2].spin_directions[:, 1], [0.2640000, 0.4870000, 2.0390000] / norm([0.2640000, 0.4870000, 2.0390000]), atol=1e-4)
    @test isapprox(spinconfigs[2].local_magfield[:, 1], [-3.11130e-02, 8.43390e-04, 3.87480e-03], atol=1e-4)

    @test isapprox(spinconfigs[2].magmom_size[16], norm([0.0150000, 0.0200000, 0.2470000]), atol=1e-4)
    @test isapprox(spinconfigs[2].spin_directions[:, 16], [0.0150000, 0.0200000, 0.2470000] / norm([0.0150000, 0.0200000, 0.2470000]), atol=1e-4)
    @test isapprox(spinconfigs[2].local_magfield[:, 16], [4.59820e-03, -6.10690e-02, 4.73170e-03], atol=1e-4)

    @test isapprox(spinconfigs[40].magmom_size[1], norm([-0.2760000, 0.3050000, 2.0750000]), atol=1e-4)
    @test isapprox(spinconfigs[40].spin_directions[:, 1], [-0.2760000, 0.3050000, 2.0750000] / norm([-0.2760000, 0.3050000, 2.0750000]), atol=1e-4)
    @test isapprox(spinconfigs[40].local_magfield[:, 1], [3.45100e-02, -1.48800e-02, 6.85560e-03], atol=1e-4)

    @test isapprox(spinconfigs[40].magmom_size[16], norm([-0.0150000, -0.0340000, 0.2390000]), atol=1e-4)
    @test isapprox(spinconfigs[40].spin_directions[:, 16], [-0.0150000, -0.0340000, 0.2390000] / norm([-0.0150000, -0.0340000, 0.2390000]), atol=1e-4)
    @test isapprox(spinconfigs[40].local_magfield[:, 16], [8.67710e-02, 1.21410e-01, 2.38430e-02], atol=1e-4)

    # @show spinconfigs[1]
end
