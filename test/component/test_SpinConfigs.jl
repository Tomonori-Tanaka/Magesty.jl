using LinearAlgebra

@testset "SpinConfig" begin
    # Test reading EMBSET file
    embset_path = normpath(joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2", "EMBSET"))
    spinconfigs::Vector{SpinConfig} = read_embset(embset_path)

    # Test first configuration
    @test length(spinconfigs) >= 1  # At least 1 configurations in the file
    @test isapprox(spinconfigs[1].energy, -121.37424)
    @test isapprox(spinconfigs[2].energy, -121.37488)

    # Test first atom in first configuration
    @test isapprox(spinconfigs[1].magmom_size[1], norm([0.1870000, 0.0000000, 2.1670000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 1], [0.1870000, 0.0000000, 2.1670000] / norm([0.1870000, 0.0000000, 2.1670000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 1], [-1.17370e-02, 4.72560e-07, 1.01990e-03], atol=1e-4)

    # Test second atom in first configuration
    @test isapprox(spinconfigs[1].magmom_size[2], norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 2], [0.0000000, 0.0000000, 2.1760000] / norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 2], [-7.03430e-04, -5.49880e-06, -2.47420e-08], atol=1e-4)

    # Test ninth atom in first configuration (Pt atom)
    @test isapprox(spinconfigs[1].magmom_size[9], norm([0.0000000, 0.0000000, 0.2380000]), atol=1e-4)
    @test isapprox(spinconfigs[1].spin_directions[:, 9], [0.0000000, 0.0000000, 0.2380000] / norm([0.0000000, 0.0000000, 0.2380000]), atol=1e-4)
    @test isapprox(spinconfigs[1].local_magfield[:, 9], [9.35110e-03, 8.61050e-05, -4.37250e-06], atol=1e-4)

    # Test first atom in second configuration
    @test isapprox(spinconfigs[2].magmom_size[1], norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
    @test isapprox(spinconfigs[2].spin_directions[:, 1], [0.0000000, 0.0000000, 2.1760000] / norm([0.0000000, 0.0000000, 2.1760000]), atol=1e-4)
    @test isapprox(spinconfigs[2].local_magfield[:, 1], [1.01860e-03, -6.52810e-06, -5.18810e-08], atol=1e-4)

    # Test local_magfield_vertical calculation
    for config in spinconfigs
        for i in 1:size(config.local_magfield, 2)
            # Check that local_magfield_vertical is perpendicular to spin_directions
            @test isapprox(dot(config.local_magfield_vertical[:, i], config.spin_directions[:, i]), 0.0, atol=1e-10)

            # Check that local_magfield_vertical is correctly calculated
            proj = dot(config.local_magfield[:, i], config.spin_directions[:, i]) * config.spin_directions[:, i]
            expected_vertical = config.local_magfield[:, i] - proj
            @test isapprox(config.local_magfield_vertical[:, i], expected_vertical, atol=1e-10)
        end
    end

    @testset "constructor input validation" begin
        good_magmom_magnitude = [1.0, 1.0]
        good_field = zeros(3, 2)

        # spin_directions must be 3 x num_atoms (row count == 3).
        @test_throws ArgumentError SpinConfig(0.0, good_magmom_magnitude, zeros(2, 2), good_field)
        @test_throws ArgumentError SpinConfig(0.0, good_magmom_magnitude, zeros(4, 2), good_field)

        # local_magfield must be 3 x num_atoms.
        @test_throws ArgumentError SpinConfig(0.0, good_magmom_magnitude, zeros(3, 2), zeros(2, 2))
        @test_throws ArgumentError SpinConfig(0.0, good_magmom_magnitude, zeros(3, 2), zeros(4, 2))

        # Column count must match length(magmom_size).
        @test_throws ArgumentError SpinConfig(0.0, good_magmom_magnitude, zeros(3, 3), good_field)

        # Negative magnetic moment sizes are rejected.
        @test_throws ArgumentError SpinConfig(0.0, [1.0, -0.5], zeros(3, 2), good_field)
    end

    @testset "spin_directions unit-norm validation" begin
        magmom_magnitude = [1.0, 1.0]
        field = zeros(3, 2)

        # A non-unit column (norm = 2.0) is rejected at the default tolerance.
        non_unit_directions = hcat([1.0, 0.0, 0.0], [2.0, 0.0, 0.0])
        @test_throws ArgumentError SpinConfig(0.0, magmom_magnitude, non_unit_directions, field)

        # A NaN column (e.g. a zero-moment row normalized as moment / 0) is
        # rejected explicitly rather than propagating NaN through predictions.
        nan_directions = hcat([1.0, 0.0, 0.0], [NaN, NaN, NaN])
        @test_throws ArgumentError SpinConfig(0.0, magmom_magnitude, nan_directions, field)

        # A column within 5e-7 of unit norm is accepted under the default 1e-6 tol.
        near_unit_offset = 5e-7
        near_unit_directions = hcat([1.0 + near_unit_offset, 0.0, 0.0], [0.0, 1.0, 0.0])
        @test SpinConfig(0.0, magmom_magnitude, near_unit_directions, field) isa SpinConfig

        # A loose tolerance accepts norms that the default would reject.
        loose_directions = hcat([1.3, 0.0, 0.0], [0.0, 1.0, 0.0])
        @test_throws ArgumentError SpinConfig(0.0, magmom_magnitude, loose_directions, field)
        @test SpinConfig(
            0.0, magmom_magnitude, loose_directions, field;
            atol_unit_norm = 0.5,
        ) isa SpinConfig
    end

    @testset "zero-moment atom in EMBSET" begin
        # A non-magnetic site reports a zero moment. `read_embset` must not
        # normalize 0/0 (which would yield a NaN direction and abort the read);
        # it stores a placeholder unit direction, and the zero magnitude makes
        # the torque vanish.
        embset = tempname()
        open(embset, "w") do io
            println(io, "# one config, atom 1 non-magnetic")
            println(io, "-1.5")
            println(io, "1  0.0 0.0 0.0   0.01 0.0 0.0")
            println(io, "2  0.0 0.0 1.5   0.02 0.0 0.0")
        end
        configs = read_embset(embset)
        rm(embset)

        @test length(configs) == 1
        sc = configs[1]
        # Zero-moment atom: magnitude is zero, direction is a unit placeholder,
        # torque is exactly zero.
        @test sc.magmom_size[1] == 0.0
        @test isapprox(norm(sc.spin_directions[:, 1]), 1.0; atol = 1e-12)
        @test all(sc.torques[:, 1] .== 0.0)
        # Magnetic atom is read normally.
        @test isapprox(sc.magmom_size[2], 1.5; atol = 1e-10)
        @test isapprox(sc.spin_directions[:, 2], [0.0, 0.0, 1.0]; atol = 1e-10)
    end

end
