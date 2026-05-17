using Test
using Magesty
using AtomsBase
using Unitful
using OffsetArrays

# Direct access to the adapter internals. Only `system_to_specs` and
# `kwargs_to_specs` are exported; the helpers `_to_angstrom`,
# `_kind_name`, `_kind_tables`, `_lookup_species`, `_lookup_pair`, and
# `_warn_if_odd` are private, but each carries a small invariant that
# the public path relies on, so each gets a direct check.
const _Adapter = Magesty.AtomsBaseAdapter

# Build a minimal periodic system with the given per-atom kinds. The
# helper keeps every test deterministic about how the AtomsBase system
# is laid out so the expected values below can be derived independently.
function _make_system(
    elements::AbstractVector{Symbol},
    cart_positions::AbstractVector;
    atom_names = nothing,
    lattice_diag::Real = 4.0,
)
    n = length(elements)
    @assert length(cart_positions) == n
    if atom_names !== nothing
        @assert length(atom_names) == n
    end
    atoms = map(1:n) do i
        pos = cart_positions[i] .* u"Å"
        if atom_names === nothing
            return Atom(elements[i], pos)
        else
            return Atom(elements[i], pos; atom_name = atom_names[i])
        end
    end
    cv = (
        [lattice_diag, 0.0, 0.0]u"Å",
        [0.0, lattice_diag, 0.0]u"Å",
        [0.0, 0.0, lattice_diag]u"Å",
    )
    return periodic_system(collect(atoms), cv)
end

@testset "AtomsBaseAdapter" begin
    # ---- _to_angstrom: unit normalization -------------------------------

    @testset "_to_angstrom" begin
        # Bare reals (no units) are taken to already be in angstrom. The
        # cutoff sentinel `-1.0` must pass through unchanged.
        @test _Adapter._to_angstrom(1.0) === 1.0
        @test _Adapter._to_angstrom(-1.0) === -1.0
        @test _Adapter._to_angstrom(0) === 0.0  # Int promotes to Float64

        # Unitful lengths convert into angstrom. All checks below use
        # decimal-prefix units whose Å values are exact in Float64.
        @test _Adapter._to_angstrom(1.0u"Å") === 1.0
        @test _Adapter._to_angstrom(1.0u"nm") === 10.0   # 1 nm  = 10 Å
        @test _Adapter._to_angstrom(1.0u"pm") === 0.01   # 1 pm  = 0.01 Å
        @test _Adapter._to_angstrom(2.5u"Å") === 2.5

        # Return type is Float64 for both branches.
        @test _Adapter._to_angstrom(1.0u"nm") isa Float64
        @test _Adapter._to_angstrom(3) isa Float64
    end

    # ---- _kind_name: per-atom kind extraction ---------------------------

    @testset "_kind_name" begin
        sys = _make_system([:Fe, :Pt], [[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]])
        # Without `atom_name`, the element symbol is the kind name.
        @test _Adapter._kind_name(sys, 1) === :Fe
        @test _Adapter._kind_name(sys, 2) === :Pt

        sys2 = _make_system(
            [:Fe, :Fe], [[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]];
            atom_names = (:Fe_4a, :Fe_8e),
        )
        # `atom_name` overrides the element symbol.
        @test _Adapter._kind_name(sys2, 1) === :Fe_4a
        @test _Adapter._kind_name(sys2, 2) === :Fe_8e
    end

    # ---- _kind_tables: kd_name, kd_int_list, kind_to_element ------------

    @testset "_kind_tables" begin
        sys = _make_system(
            [:Fe, :Pt, :Fe], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
        )
        kd_name, kd_int_list, kind_to_element = _Adapter._kind_tables(sys)
        # Kind order follows first appearance.
        @test kd_name == ["Fe", "Pt"]
        @test kd_int_list == [1, 2, 1]
        @test kind_to_element == Dict(:Fe => :Fe, :Pt => :Pt)

        sys_sub = _make_system(
            [:Fe, :Fe, :Pt],
            [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
            atom_names = (:Fe_4a, :Fe_8e, :Pt_8f),
        )
        kd_name_sub, kd_int_list_sub, kind_to_element_sub =
            _Adapter._kind_tables(sys_sub)
        # Sublabels distinct from the element symbol, and from each other.
        @test kd_name_sub == ["Fe_4a", "Fe_8e", "Pt_8f"]
        @test kd_int_list_sub == [1, 2, 3]
        @test kind_to_element_sub == Dict(
            :Fe_4a => :Fe, :Fe_8e => :Fe, :Pt_8f => :Pt,
        )
    end

    # ---- _lookup_species: kind preferred, element fallback --------------

    @testset "_lookup_species" begin
        # Kind key present: take the kind value, ignore element fallback.
        d_kind = Dict(:Fe_4a => 3, :Fe => 2)
        @test _Adapter._lookup_species(d_kind, :Fe_4a, :Fe, "lmax") == 3

        # Kind key absent, element key present: element fan-out.
        d_elem = Dict(:Fe => 2)
        @test _Adapter._lookup_species(d_elem, :Fe_4a, :Fe, "lmax") == 2

        # Neither key present: error message names the kind and element.
        d_empty = Dict{Symbol, Int}(:Co => 4)
        err = ArgumentError("lmax not specified for kind :Fe_4a (element :Fe)")
        @test_throws err _Adapter._lookup_species(d_empty, :Fe_4a, :Fe, "lmax")
    end

    # ---- _lookup_pair: both orderings, element fan-out ------------------

    @testset "_lookup_pair" begin
        # Kind pair, original ordering present.
        d_kind_fwd = Dict((:Fe, :Pt) => 5.0)
        @test _Adapter._lookup_pair(d_kind_fwd, :Fe, :Pt, :Fe, :Pt) == 5.0
        @test _Adapter._lookup_pair(d_kind_fwd, :Pt, :Fe, :Pt, :Fe) == 5.0  # reverse

        # Kind absent; element pair (forward) present.
        d_elem_fwd = Dict((:Fe, :Pt) => 5.0)
        @test _Adapter._lookup_pair(
            d_elem_fwd, :Fe_4a, :Pt_8f, :Fe, :Pt,
        ) == 5.0

        # Kind absent; element pair (reverse) present.
        d_elem_rev = Dict((:Pt, :Fe) => 4.0)
        @test _Adapter._lookup_pair(
            d_elem_rev, :Fe_4a, :Pt_8f, :Fe, :Pt,
        ) == 4.0

        # Nothing matches.
        d_empty = Dict{Tuple{Symbol, Symbol}, Float64}()
        @test_throws ArgumentError _Adapter._lookup_pair(
            d_empty, :Fe, :Pt, :Fe, :Pt,
        )
    end

    # ---- _warn_if_odd ---------------------------------------------------

    @testset "_warn_if_odd" begin
        # Even values pass through silently and return unchanged.
        @test _Adapter._warn_if_odd(0, "lmax") == 0
        @test _Adapter._warn_if_odd(2, "lmax") == 2
        @test _Adapter._warn_if_odd(4, "lmax") == 4

        # Odd values emit a `@warn` and still return the input.
        @test_logs (:warn, r"lmax.*= 3 is odd") begin
            @test _Adapter._warn_if_odd(3, "lmax") == 3
        end
    end

    # ---- system_to_specs: end-to-end AtomsBase path ---------------------

    @testset "system_to_specs end-to-end" begin
        # 2 Fe atoms at (0,0,0) and (1.5,1.5,1.5) Å in a 3 Å cubic cell.
        # Expected fractional positions: (0,0,0) and (0.5,0.5,0.5).
        sys = _make_system(
            [:Fe, :Fe], [[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]];
            lattice_diag = 3.0,
        )
        interaction = (
            body1 = (lmax = Dict(:Fe => 2),),
            body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0)),
        )
        system, inter, options = _Adapter.system_to_specs(
            sys, interaction; name = "test_system",
        )

        # SystemSpec fields.
        @test system.name == "test_system"
        @test system.num_atoms == 2
        @test system.kd_name == ["Fe"]
        @test system.kd_int_list == [1, 1]
        @test system.lattice_vectors == [
            3.0 0.0 0.0;
            0.0 3.0 0.0;
            0.0 0.0 3.0
        ]
        @test system.x_fractional ≈ [
            0.0 0.5;
            0.0 0.5;
            0.0 0.5
        ] atol = 1e-12
        @test system.is_periodic == [true, true, true]

        # InteractionSpec fields. The cutoff sentinel `-1.0` passes
        # through unchanged.
        @test inter.nbody == 2
        @test inter.body1_lmax == [2]
        @test inter.bodyn_lsum[2] == 2
        @test inter.bodyn_cutoff[2, 1, 1] == -1.0

        # SymmetryOptions defaults (tolerance_sym = 1e-3, isotropy = false).
        @test options.tolerance_sym == 1e-3
        @test options.isotropy == false
    end

    @testset "system_to_specs respects keyword overrides" begin
        sys = _make_system([:Fe], [[0.0, 0.0, 0.0]]; lattice_diag = 3.0)
        interaction = (
            body1 = (lmax = Dict(:Fe => 0),),
            body2 = (lsum = 0, cutoff = Dict((:Fe, :Fe) => -1.0)),
        )
        _, _, options = _Adapter.system_to_specs(
            sys, interaction;
            name = "custom",
            tolerance_sym = 1e-5,
            isotropy = true,
        )
        @test options.tolerance_sym == 1e-5
        @test options.isotropy == true
    end

    @testset "system_to_specs handles Unitful cutoff" begin
        sys = _make_system([:Fe, :Fe], [[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]])
        interaction = (
            body1 = (lmax = Dict(:Fe => 0),),
            # 1 nm cutoff -> 10 Å when normalized.
            body2 = (lsum = 0, cutoff = Dict((:Fe, :Fe) => 1.0u"nm")),
        )
        _, inter, _ = _Adapter.system_to_specs(sys, interaction)
        @test inter.bodyn_cutoff[2, 1, 1] ≈ 10.0
    end

    @testset "system_to_specs element fan-out across sublabels" begin
        sys = _make_system(
            [:Fe, :Fe], [[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]];
            atom_names = (:Fe_4a, :Fe_8e),
        )
        # An element-only :Fe key must fan out to both Fe_4a and Fe_8e.
        interaction = (
            body1 = (lmax = Dict(:Fe => 2),),
            body2 = (lsum = 2, cutoff = Dict((:Fe, :Fe) => -1.0)),
        )
        system, inter, _ = _Adapter.system_to_specs(sys, interaction)
        @test Set(system.kd_name) == Set(["Fe_4a", "Fe_8e"])
        # body1_lmax has one entry per kind; both fan out to the element value 2.
        @test length(inter.body1_lmax) == 2
        @test all(inter.body1_lmax .== 2)
        # The 2x2 cutoff matrix is filled symmetrically with -1.0.
        @test all(inter.bodyn_cutoff[2, i, j] == -1.0 for i in 1:2, j in 1:2)
    end

    @testset "system_to_specs odd lsum emits a warning" begin
        sys = _make_system([:Fe], [[0.0, 0.0, 0.0]])
        interaction = (
            body1 = (lmax = Dict(:Fe => 0),),
            body2 = (lsum = 3, cutoff = Dict((:Fe, :Fe) => -1.0)),
        )
        @test_logs (:warn, r"lsum for body2 = 3 is odd") begin
            _Adapter.system_to_specs(sys, interaction)
        end
    end

    # ---- kwargs_to_specs: keyword path (no AtomsBase) -------------------

    @testset "kwargs_to_specs end-to-end" begin
        system, inter, options = _Adapter.kwargs_to_specs(
            lattice = [3.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 3.0],
            kd = [:X],
            kd_list = [1, 1],
            positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            periodicity = (true, true, true),
            interaction = (
                body1 = (lmax = Dict(:X => 0),),
                body2 = (lsum = 2, cutoff = Dict((:X, :X) => -1.0)),
            ),
            name = "kw",
        )
        @test system.name == "kw"
        @test system.num_atoms == 2
        @test system.kd_name == ["X"]
        @test system.kd_int_list == [1, 1]
        @test system.lattice_vectors == [
            3.0 0.0 0.0;
            0.0 3.0 0.0;
            0.0 0.0 3.0
        ]
        # Keyword positions are already fractional — no `inv(lattice) * cart`
        # conversion.
        @test system.x_fractional == [
            0.0 0.5;
            0.0 0.5;
            0.0 0.5
        ]
        @test system.is_periodic == [true, true, true]
        @test inter.nbody == 2
        @test inter.body1_lmax == [0]
        @test inter.bodyn_lsum[2] == 2
        @test inter.bodyn_cutoff[2, 1, 1] == -1.0
        @test options.tolerance_sym == 1e-3
        @test options.isotropy == false
    end

    @testset "kwargs_to_specs validates lattice shape" begin
        # 2x3 lattice → ArgumentError.
        @test_throws ArgumentError _Adapter.kwargs_to_specs(
            lattice = [1.0 0.0 0.0; 0.0 1.0 0.0],
            kd = [:X],
            kd_list = [1],
            positions = [[0.0, 0.0, 0.0]],
            interaction = (
                body1 = (lmax = Dict(:X => 0),),
            ),
        )
    end

    @testset "kwargs_to_specs accepts Unitful lattice entries" begin
        # 1 nm lattice in each direction.
        system, _, _ = _Adapter.kwargs_to_specs(
            lattice = [1.0u"nm" 0.0u"nm" 0.0u"nm";
                       0.0u"nm" 1.0u"nm" 0.0u"nm";
                       0.0u"nm" 0.0u"nm" 1.0u"nm"],
            kd = [:X],
            kd_list = [1],
            positions = [[0.0, 0.0, 0.0]],
            interaction = (
                body1 = (lmax = Dict(:X => 0),),
            ),
        )
        @test system.lattice_vectors == [
            10.0 0.0 0.0;
            0.0 10.0 0.0;
            0.0 0.0 10.0
        ]
    end

    @testset "kwargs_to_specs rejects missing body1" begin
        @test_throws ArgumentError _Adapter.kwargs_to_specs(
            lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
            kd = [:X],
            kd_list = [1],
            positions = [[0.0, 0.0, 0.0]],
            # No body1 in interaction -> _build_interaction_spec throws.
            interaction = (body2 = (lsum = 2, cutoff = Dict((:X, :X) => -1.0)),),
        )
    end

    @testset "kwargs_to_specs rejects missing bodyN when nbody >= N" begin
        # nbody is inferred from `length(interaction)`. Passing body1 + body3
        # (no body2) -> nbody = 2 but `body2` key missing -> ArgumentError.
        @test_throws ArgumentError _Adapter.kwargs_to_specs(
            lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
            kd = [:X],
            kd_list = [1],
            positions = [[0.0, 0.0, 0.0]],
            interaction = (
                body1 = (lmax = Dict(:X => 0),),
                body3 = (lsum = 2, cutoff = Dict((:X, :X) => -1.0)),
            ),
        )
    end
end
