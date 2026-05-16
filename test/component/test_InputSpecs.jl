using .InputSpecs
using OffsetArrays
using Test

@testset "InputSpecs" begin
    # Reusable fixtures.
    valid_lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    valid_positions = [0.0 0.5; 0.0 0.5; 0.0 0.5]
    valid_kd = ["H", "O"]
    valid_kd_int = [1, 2]

    @testset "SystemSpec validation" begin
        # Happy path.
        sys = SystemSpec(
            name = "water",
            num_atoms = 2,
            kd_name = valid_kd,
            kd_int_list = valid_kd_int,
            lattice_vectors = valid_lattice,
            x_fractional = valid_positions,
        )
        @test sys.name == "water"
        @test sys.num_atoms == 2
        @test sys.kd_name == valid_kd
        @test sys.kd_int_list == valid_kd_int
        @test sys.lattice_vectors == valid_lattice
        @test sys.x_fractional == valid_positions
        @test sys.is_periodic == [true, true, true]   # default

        # Reject empty name.
        @test_throws ArgumentError SystemSpec(
            name = "",
            num_atoms = 2, kd_name = valid_kd, kd_int_list = valid_kd_int,
            lattice_vectors = valid_lattice, x_fractional = valid_positions,
        )

        # Reject non-positive atom count.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 0,
            kd_name = valid_kd, kd_int_list = Int[],
            lattice_vectors = valid_lattice,
            x_fractional = zeros(3, 0),
        )

        # Reject empty species list.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 1,
            kd_name = String[], kd_int_list = [1],
            lattice_vectors = valid_lattice,
            x_fractional = zeros(3, 1),
        )

        # Reject duplicate species.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 1,
            kd_name = ["H", "H"], kd_int_list = [1],
            lattice_vectors = valid_lattice,
            x_fractional = zeros(3, 1),
        )

        # Reject non-3x3 lattice.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 2,
            kd_name = valid_kd, kd_int_list = valid_kd_int,
            lattice_vectors = [1.0 0.0; 0.0 1.0],
            x_fractional = valid_positions,
        )

        # Reject mismatched x_fractional column count.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 2,
            kd_name = valid_kd, kd_int_list = valid_kd_int,
            lattice_vectors = valid_lattice,
            x_fractional = zeros(3, 1),
        )

        # Reject mismatched kd_int_list length.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 2,
            kd_name = valid_kd, kd_int_list = [1],
            lattice_vectors = valid_lattice,
            x_fractional = valid_positions,
        )

        # Reject kd_int out of range.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 2,
            kd_name = valid_kd, kd_int_list = [1, 5],
            lattice_vectors = valid_lattice,
            x_fractional = valid_positions,
        )

        # Reject is_periodic of wrong length.
        @test_throws ArgumentError SystemSpec(
            name = "x", num_atoms = 2,
            kd_name = valid_kd, kd_int_list = valid_kd_int,
            lattice_vectors = valid_lattice,
            x_fractional = valid_positions,
            is_periodic = [true, true],
        )
    end

    @testset "InteractionSpec validation, nbody=1" begin
        spec = InteractionSpec(
            nbody = 1,
            body1_lmax = [2, 3],
            bodyn_lsum = OffsetArray(Int[], 2:1),
            bodyn_cutoff = OffsetArray(zeros(Float64, 0, 2, 2), 2:1, 1:2, 1:2),
            kd_name = valid_kd,
        )
        @test spec.nbody == 1
        @test spec.body1_lmax == [2, 3]
        @test length(spec.bodyn_lsum) == 0
        @test size(parent(spec.bodyn_cutoff)) == (0, 2, 2)

        # Reject non-empty bodyn_lsum when nbody == 1.
        @test_throws ArgumentError InteractionSpec(
            nbody = 1, body1_lmax = [2, 3],
            bodyn_lsum = OffsetArray([4], 2:2),
            bodyn_cutoff = OffsetArray(zeros(Float64, 0, 2, 2), 2:1, 1:2, 1:2),
            kd_name = valid_kd,
        )

        # Reject body1_lmax length mismatch.
        @test_throws ArgumentError InteractionSpec(
            nbody = 1, body1_lmax = [2],
            bodyn_lsum = OffsetArray(Int[], 2:1),
            bodyn_cutoff = OffsetArray(zeros(Float64, 0, 2, 2), 2:1, 1:2, 1:2),
            kd_name = valid_kd,
        )
    end

    @testset "InteractionSpec validation, nbody>=2" begin
        nkd = 2
        cutoff = OffsetArray(zeros(Float64, 1, nkd, nkd), 2:2, 1:nkd, 1:nkd)
        cutoff[2, 1, 1] = 2.0
        cutoff[2, 1, 2] = 3.0
        cutoff[2, 2, 1] = 3.0
        cutoff[2, 2, 2] = 4.0
        lsum = OffsetArray([4], 2:2)

        spec = InteractionSpec(
            nbody = 2,
            body1_lmax = [1, 3],
            bodyn_lsum = lsum,
            bodyn_cutoff = cutoff,
            kd_name = valid_kd,
        )
        @test spec.bodyn_lsum[2] == 4
        @test spec.bodyn_cutoff[2, 1, 2] ≈ 3.0
        @test spec.bodyn_cutoff[2, 2, 1] ≈ 3.0

        # Reject asymmetric cutoff.
        asym = OffsetArray(zeros(Float64, 1, nkd, nkd), 2:2, 1:nkd, 1:nkd)
        asym[2, 1, 1] = 2.0
        asym[2, 1, 2] = 3.0
        asym[2, 2, 1] = 9.9                  # break symmetry
        asym[2, 2, 2] = 4.0
        @test_throws ArgumentError InteractionSpec(
            nbody = 2, body1_lmax = [1, 3],
            bodyn_lsum = lsum, bodyn_cutoff = asym,
            kd_name = valid_kd,
        )

        # Reject nbody = 0.
        @test_throws ArgumentError InteractionSpec(
            nbody = 0, body1_lmax = Int[],
            bodyn_lsum = OffsetArray(Int[], 2:1),
            bodyn_cutoff = OffsetArray(zeros(Float64, 0, 0, 0), 2:1, 1:0, 1:0),
            kd_name = String[],
        )

        # Reject wrong bodyn_cutoff axes.
        wrong_axes = OffsetArray(zeros(Float64, 2, nkd, nkd), 2:3, 1:nkd, 1:nkd)
        @test_throws ArgumentError InteractionSpec(
            nbody = 2, body1_lmax = [1, 3],
            bodyn_lsum = lsum, bodyn_cutoff = wrong_axes,
            kd_name = valid_kd,
        )
    end

    @testset "SymmetryOptions defaults and validation" begin
        opts = SymmetryOptions()
        @test opts.tolerance_sym ≈ 1e-3
        @test opts.isotropy == false

        opts2 = SymmetryOptions(tolerance_sym = 1e-5, isotropy = true)
        @test opts2.tolerance_sym ≈ 1e-5
        @test opts2.isotropy == true

        @test_throws ArgumentError SymmetryOptions(tolerance_sym = 0.0)
        @test_throws ArgumentError SymmetryOptions(tolerance_sym = -1.0)
    end

    @testset "expand_pair_table specificity resolution" begin
        kd3 = ["Fe", "Ni", "Co"]

        # "*-*" alone fills every pair.
        m = InputSpecs.expand_pair_table(kd3, Dict("*-*" => 8.0))
        @test all(m .== 8.0)
        @test size(m) == (3, 3)
        for i in 1:3, j in 1:3
            @test m[i, j] == m[j, i]
        end

        # "*-*" + "Fe-Fe": Fe-Fe wins for that pair, others stay 8.0.
        m = InputSpecs.expand_pair_table(kd3, Dict("*-*" => 8.0, "Fe-Fe" => 12.0))
        fe = findfirst(==("Fe"), kd3)
        @test m[fe, fe] ≈ 12.0
        for i in 1:3, j in 1:3
            (i == fe && j == fe) && continue
            @test m[i, j] ≈ 8.0
        end

        # "*-*" + "Fe-*": every pair touching Fe goes to 10.0, others 8.0.
        m = InputSpecs.expand_pair_table(kd3, Dict("*-*" => 8.0, "Fe-*" => 10.0))
        for i in 1:3, j in 1:3
            if i == fe || j == fe
                @test m[i, j] ≈ 10.0
            else
                @test m[i, j] ≈ 8.0
            end
        end

        # All-concrete schema (legacy, no wildcards): every required pair listed.
        m = InputSpecs.expand_pair_table(
            ["H", "O"],
            Dict("H-H" => 2.0, "H-O" => 3.0, "O-O" => 4.0),
        )
        @test m[1, 1] ≈ 2.0
        @test m[1, 2] ≈ 3.0
        @test m[2, 1] ≈ 3.0
        @test m[2, 2] ≈ 4.0

        # Unordered pair: "Ni-Fe" equivalent to "Fe-Ni" — supplying only one is enough.
        m = InputSpecs.expand_pair_table(
            ["Fe", "Ni"],
            Dict("Fe-Fe" => 1.0, "Ni-Fe" => 2.0, "Ni-Ni" => 3.0),
        )
        @test m[1, 2] ≈ 2.0
        @test m[2, 1] ≈ 2.0
    end

    @testset "expand_pair_table error conditions" begin
        kd3 = ["Fe", "Ni", "Co"]

        # Tier-1 ambiguity: "Fe-*" and "*-Ni" both cover Fe-Ni with different values.
        @test_throws ArgumentError InputSpecs.expand_pair_table(
            kd3,
            Dict("*-*" => 8.0, "Fe-*" => 10.0, "*-Ni" => 11.0),
        )

        # Duplicate equivalent keys: "Fe-Ni" and "Ni-Fe" specify the same pair.
        @test_throws ArgumentError InputSpecs.expand_pair_table(
            ["Fe", "Ni"],
            Dict("Fe-Fe" => 1.0, "Fe-Ni" => 2.0, "Ni-Fe" => 2.5, "Ni-Ni" => 3.0),
        )

        # Unknown species in a concrete key.
        @test_throws ArgumentError InputSpecs.expand_pair_table(
            kd3,
            Dict("Fe-Xx" => 10.0),
        )

        # Missing coverage: no key matches Ni-Co.
        @test_throws ArgumentError InputSpecs.expand_pair_table(
            kd3,
            Dict("Fe-Fe" => 1.0, "Fe-Ni" => 2.0, "Fe-Co" => 3.0, "Ni-Ni" => 4.0),
        )

        # Malformed key.
        @test_throws ArgumentError InputSpecs.expand_pair_table(
            ["Fe"],
            Dict("Fe" => 1.0),
        )
    end

    @testset "expand_species_table" begin
        kd3 = ["Fe", "Ni", "Co"]

        # Wildcard alone.
        v = InputSpecs.expand_species_table(kd3, Dict("*" => 2))
        @test v == [2, 2, 2]

        # Concrete overrides wildcard.
        v = InputSpecs.expand_species_table(kd3, Dict("*" => 2, "Fe" => 4))
        @test v[1] == 4
        @test v[2] == 2
        @test v[3] == 2

        # All concrete.
        v = InputSpecs.expand_species_table(kd3, Dict("Fe" => 1, "Ni" => 2, "Co" => 3))
        @test v == [1, 2, 3]

        # Empty: returns zeros (preserves legacy "body1 omitted" default).
        v = InputSpecs.expand_species_table(kd3, Dict{String, Int}())
        @test v == [0, 0, 0]

        # Missing coverage without wildcard → error.
        @test_throws ArgumentError InputSpecs.expand_species_table(
            kd3, Dict("Fe" => 1, "Ni" => 2),
        )

        # Unknown species.
        @test_throws ArgumentError InputSpecs.expand_species_table(
            kd3, Dict("Xx" => 1),
        )
    end

    @testset "parse_toml_inputs schema mapping" begin
        # The legacy TOML schema must produce specific values. Tests pin the
        # mapping from each TOML key to the corresponding typed-value field.
        input_dict = Dict{String, Any}(
            "general" => Dict{String, Any}(
                "name" => "test_system",
                "nat" => 2,
                "kd" => ["H", "O"],
                "periodicity" => [true, true, true],
            ),
            "symmetry" => Dict{String, Any}(
                "tolerance" => 1e-3,
            ),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict{String, Any}(
                    "lmax" => Dict{String, Any}("H" => 1, "O" => 3),
                ),
                "body2" => Dict{String, Any}(
                    "lsum" => 4,
                    "cutoff" => Dict{String, Any}(
                        "H-H" => 2.0,
                        "H-O" => 3.0,
                        "O-O" => 4.0,
                    ),
                ),
            ),
            "structure" => Dict{String, Any}(
                "lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "kd_list" => [1, 2],
                "position" => [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            ),
        )

        sys, inter, opts = parse_toml_inputs(input_dict)

        @test sys.name == "test_system"
        @test sys.num_atoms == 2
        @test sys.kd_name == ["H", "O"]
        @test sys.kd_int_list == [1, 2]
        @test sys.lattice_vectors == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test sys.x_fractional[:, 1] ≈ [0.0, 0.0, 0.0]
        @test sys.x_fractional[:, 2] ≈ [0.5, 0.5, 0.5]
        @test sys.is_periodic == [true, true, true]

        @test inter.nbody == 2
        @test inter.body1_lmax == [1, 3]
        @test inter.bodyn_lsum[2] == 4
        @test inter.bodyn_cutoff[2, 1, 1] ≈ 2.0
        @test inter.bodyn_cutoff[2, 1, 2] ≈ 3.0
        @test inter.bodyn_cutoff[2, 2, 1] ≈ 3.0
        @test inter.bodyn_cutoff[2, 2, 2] ≈ 4.0

        @test opts.tolerance_sym ≈ 1e-3
        @test opts.isotropy == false
    end

    @testset "parse_toml_inputs nbody=1 edge case" begin
        input_dict = Dict{String, Any}(
            "general" => Dict{String, Any}(
                "name" => "nbody1",
                "nat" => 1,
                "kd" => ["Fe"],
            ),
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}(
                "nbody" => 1,
                "body1" => Dict{String, Any}(
                    "lmax" => Dict{String, Any}("Fe" => 2),
                ),
            ),
            "structure" => Dict{String, Any}(
                "lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "kd_list" => [1],
                "position" => [[0.0, 0.0, 0.0]],
            ),
        )
        sys, inter, opts = parse_toml_inputs(input_dict)
        @test inter.nbody == 1
        @test inter.body1_lmax == [2]
        @test length(inter.bodyn_lsum) == 0
        @test size(parent(inter.bodyn_cutoff)) == (0, 1, 1)
        @test opts.tolerance_sym ≈ 1e-3   # default
    end

    @testset "parse_toml_inputs wildcards equivalent to explicit form" begin
        # Two TOML-shape inputs that should resolve to the same InteractionSpec.
        base_general = Dict{String, Any}(
            "name" => "x", "nat" => 1, "kd" => ["Fe", "Ni", "Co"],
        )
        base_structure = Dict{String, Any}(
            "lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "kd_list" => [1],
            "position" => [[0.0, 0.0, 0.0]],
        )

        explicit = Dict{String, Any}(
            "general" => base_general,
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict("lmax" => Dict("Fe" => 2, "Ni" => 2, "Co" => 2)),
                "body2" => Dict(
                    "lsum" => 4,
                    "cutoff" => Dict(
                        "Fe-Fe" => 8.0, "Fe-Ni" => 8.0, "Fe-Co" => 8.0,
                        "Ni-Ni" => 8.0, "Ni-Co" => 8.0, "Co-Co" => 8.0,
                    ),
                ),
            ),
            "structure" => base_structure,
        )
        wildcard = Dict{String, Any}(
            "general" => base_general,
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict("lmax" => Dict("*" => 2)),
                "body2" => Dict(
                    "lsum" => 4,
                    "cutoff" => Dict("*-*" => 8.0),
                ),
            ),
            "structure" => base_structure,
        )

        _, inter_explicit, _ = parse_toml_inputs(explicit)
        _, inter_wildcard, _ = parse_toml_inputs(wildcard)

        @test inter_explicit.body1_lmax == inter_wildcard.body1_lmax
        @test inter_explicit.bodyn_cutoff == inter_wildcard.bodyn_cutoff

        # And against an override case: "*-*" 8.0 + "Fe-Fe" 12.0 vs the
        # equivalent fully-explicit form.
        override = Dict{String, Any}(
            "general" => base_general,
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict("lmax" => Dict("*" => 2)),
                "body2" => Dict(
                    "lsum" => 4,
                    "cutoff" => Dict("*-*" => 8.0, "Fe-Fe" => 12.0),
                ),
            ),
            "structure" => base_structure,
        )
        override_explicit = Dict{String, Any}(
            "general" => base_general,
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict("lmax" => Dict("*" => 2)),
                "body2" => Dict(
                    "lsum" => 4,
                    "cutoff" => Dict(
                        "Fe-Fe" => 12.0, "Fe-Ni" => 8.0, "Fe-Co" => 8.0,
                        "Ni-Ni" => 8.0, "Ni-Co" => 8.0, "Co-Co" => 8.0,
                    ),
                ),
            ),
            "structure" => base_structure,
        )
        _, inter_or, _ = parse_toml_inputs(override)
        _, inter_oe, _ = parse_toml_inputs(override_explicit)
        @test inter_or.bodyn_cutoff == inter_oe.bodyn_cutoff
    end

    @testset "parse_toml_inputs missing section" begin
        bad = Dict{String, Any}(
            "symmetry" => Dict{String, Any}(),
            "interaction" => Dict{String, Any}("nbody" => 1),
            "structure" => Dict{String, Any}(
                "lattice" => [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
                "kd_list" => [1], "position" => [[0.0,0.0,0.0]],
            ),
        )
        @test_throws ArgumentError parse_toml_inputs(bad)
    end
end
