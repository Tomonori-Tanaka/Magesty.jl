using Test
using Magesty

# Regression test for the `l_vec == [1, 3]` (and `[3, 1]`) workaround
# that used to live in `SALCBases.listup_coupled_basislist`.
#
# History:
# - Commit 8b389b4 (Feb 2026) added `if sort(l_vec) == [1, 3]; continue; end`
#   as an exceptional handling for a bug in the projection matrix construction.
# - Commit 3fc4180 (Mar 2026) refactored `find_translation_atoms` and
#   commented the guard out (the bug was presumed fixed). No regression
#   test was added then.
# - This file (R3 cleanup) adds the missing coverage: a 2-atom BCC cell
#   with `body2.lsum = 4` exercises `l_vec in {[1,3], [2,2], [3,1]}`
#   through the full `SALCBasis` build, including
#   `projection_matrix_coupled_basis`. If the bug ever resurfaces, this
#   test fails before the cleanup commit can be merged.
#
# See DESIGN_NOTES.md section "R3" for the design memo.

@testset "SALCBases l_vec == [1, 3] regression" begin
    @testset "listup_coupled_basislist enumerates [1,3]" begin
        cb_list = Magesty.SALCBases.listup_coupled_basislist(
            [1, 2], 4; isotropy = false,
        )
        @test !isempty(cb_list)
        # At least one composition with sort(l_vec) == [1, 3] must be present.
        l_vecs = [collect(cb.ls) for cb in cb_list]
        @test any(v -> sort(v) == [1, 3], l_vecs)
        # The [2, 2] composition (also reached at lsum=4) must be present too,
        # to confirm we are not silently dropping any compositions.
        @test any(v -> sort(v) == [2, 2], l_vecs)
    end

    @testset "Full SALCBasis build on BCC with body2.lsum = 4" begin
        # Minimal conventional BCC cell (2 atoms: corner + body center)
        # with body2.lsum = 4. The projection_matrix_coupled_basis routine
        # must complete for the [1,3] / [2,2] / [3,1] compositions without
        # erroring.
        input_dict = Dict{String, Any}(
            "general" => Dict{String, Any}(
                "name" => "bcc_lsum4",
                "kd" => ["X"],
                "nat" => 2,
                "periodicity" => [true, true, true],
            ),
            "symmetry" => Dict{String, Any}(
                "tolerance" => 1e-5,
                "isotropy" => true,
            ),
            "interaction" => Dict{String, Any}(
                "nbody" => 2,
                "body1" => Dict{String, Any}(
                    "lmax" => Dict{String, Any}("X" => 0),
                ),
                "body2" => Dict{String, Any}(
                    "cutoff" => Dict{String, Any}("X-X" => -1),
                    "lsum" => 4,
                ),
            ),
            "structure" => Dict{String, Any}(
                "kd_list" => [1, 1],
                "lattice" => [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                "position" => [
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.5],
                ],
            ),
        )

        system = Magesty.System(input_dict; verbosity = false)
        @test length(system.basisset.salc_list) > 0
    end
end
