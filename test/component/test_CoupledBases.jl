@testset "CoupledBases" begin
    using .AngularMomentumCoupling
    using .CoupledBases: CoupledBasis, CoupledBasis_with_coefficient,
        tesseral_coupled_bases_from_tesseral_bases, reorder_atoms,
        convert_to_coupled_basis
    using Test

    @testset "CoupledBasis" begin
        @testset "Constructor" begin
            atoms = [10, 20]
            ls = [1, 1]
            cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
            @test length(cb_list) == 3
            # Lf = 0
            cb_lf0 = cb_list[1]
            @test cb_lf0.ls == [1, 1]
            @test cb_lf0.Lf == 0
            @test cb_lf0.Lseq == Int[]
            @test cb_lf0.atoms == [10, 20]
            @test isapprox(
                cb_lf0.coeff_tensor,
                [1/sqrt(3) 0.0 0.0; 0.0 1/sqrt(3) 0.0; 0.0 0.0 1/sqrt(3)],
                atol = 1e-10,
            )
        end

        @testset "reorder_atoms" begin
            @testset "Two-body case" begin
                atoms = [1, 2]
                ls = [3, 1]
                cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
                @test length(cb_list) == 3

                cb = cb_list[1]
                original_atoms = copy(cb.atoms)
                original_ls = copy(cb.ls)
                original_tensor = copy(cb.coeff_tensor)
                @test size(original_tensor) == (7, 3, 5)

                # Test with new_atoms that need sorting: [2, 1] -> [1, 2]
                new_atoms = [4, 1]
                cb_new = reorder_atoms(cb, new_atoms)

                # Check that atoms are sorted
                @test cb_new.atoms == [1, 4]
                # Check that ls is permuted accordingly
                @test cb_new.ls == [1, 3]
                # Check that Lf and Lseq are unchanged
                @test cb_new.Lf == cb.Lf
                @test cb_new.Lseq == cb.Lseq
                @test size(cb_new.coeff_tensor) == (3, 7, 5)
            end

            @testset "Three-body case" begin
                atoms = [1, 2, 3]
                ls = [1, 1, 1]
                cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
                @test length(cb_list) > 0

                cb = cb_list[1]
                original_atoms = copy(cb.atoms)
                original_ls = copy(cb.ls)
                original_tensor = copy(cb.coeff_tensor)

                # Test with new_atoms = [5, 1, 2] -> sorted to [1, 2, 5]
                new_atoms = [5, 1, 2]
                cb_new = reorder_atoms(cb, new_atoms)

                # Check that atoms are sorted
                @test cb_new.atoms == [1, 2, 5]
                # Check that ls is permuted accordingly (p = [2, 3, 1])
                @test cb_new.ls == original_ls[[2, 3, 1]]
                # Check that tensor dimensions are permuted
                @test size(cb_new.coeff_tensor) == size(original_tensor)
                # Check that Lf and Lseq are unchanged
                @test cb_new.Lf == cb.Lf
                @test cb_new.Lseq == cb.Lseq
            end

            @testset "Single atom case (N=1)" begin
                atoms = [1]
                ls = [2]
                cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
                @test length(cb_list) > 0

                cb = cb_list[1]
                original_atoms = copy(cb.atoms)
                original_ls = copy(cb.ls)

                # For N=1, sorting doesn't change anything, but we can test the function
                new_atoms = [5]
                cb_new = reorder_atoms(cb, new_atoms)

                @test cb_new.atoms == [5]
                @test cb_new.ls == original_ls
                @test cb_new.Lf == cb.Lf
                @test cb_new.Lseq == cb.Lseq
                @test size(cb_new.coeff_tensor) == size(cb.coeff_tensor)
            end
        end
    end

    @testset "reorder_atoms (additional tests)" begin
        @testset "Two-body case" begin
            atoms = [1, 2]
            ls = [1, 1]
            cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
            @test length(cb_list) > 0

            cb = cb_list[1]
            original_atoms = copy(cb.atoms)
            original_ls = copy(cb.ls)
            original_tensor = copy(cb.coeff_tensor)

            # Test with new_atoms that need sorting: [5, 1, 2] -> [1, 2, 5]
            # But for two-body, we'll use [2, 1] -> [1, 2]
            new_atoms = [2, 1]
            cb_new = reorder_atoms(cb, new_atoms)

            # Check that atoms are sorted
            @test cb_new.atoms == [1, 2]
            # Check that ls is permuted accordingly
            @test cb_new.ls == original_ls[[2, 1]]
            # Check that tensor dimensions are permuted
            @test size(cb_new.coeff_tensor) == size(original_tensor)
            # Check that Lf and Lseq are unchanged
            @test cb_new.Lf == cb.Lf
            @test cb_new.Lseq == cb.Lseq

            # Lf=1
            cb = cb_list[2]
            # Test with new_atoms = [5, 1, 2] -> sorted to [1, 2, 5]
            new_atoms = [2, 1]
            cb_new = reorder_atoms(cb, new_atoms)


            function tensor_inner_product(
                tensor1::AbstractArray{T, N},
                tensor2::AbstractArray{T, N},
            ) where {T, N}
                return sum(conj.(tensor1) .* tensor2)
            end
            # @show tensor_inner_product(cb.coeff_tensor, cb_new.coeff_tensor)
        end

        @testset "Three-body case" begin
            atoms = [1, 2, 3]
            ls = [1, 1, 1]
            cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
            @test length(cb_list) > 0

            cb = cb_list[1]
            original_atoms = copy(cb.atoms)
            original_ls = copy(cb.ls)
            original_tensor = copy(cb.coeff_tensor)

            # Test with new_atoms = [5, 1, 2] -> sorted to [1, 2, 5]
            new_atoms = [5, 1, 2]
            cb_new = reorder_atoms(cb, new_atoms)

            # Check that atoms are sorted
            @test cb_new.atoms == [1, 2, 5]
            # Check that ls is permuted accordingly (p = [2, 3, 1])
            @test cb_new.ls == original_ls[[2, 3, 1]]
            # Check that tensor dimensions are permuted
            @test size(cb_new.coeff_tensor) == size(original_tensor)
            # Check that Lf and Lseq are unchanged
            @test cb_new.Lf == cb.Lf
            @test cb_new.Lseq == cb.Lseq

            # Verify tensor permutation: first 3 dimensions permuted, last dimension (Mf) unchanged
            nd = ndims(original_tensor)
            @test nd == 4  # N+1 = 3+1
            # The permutation should be [2, 3, 1, 4] for dimensions
            expected_size = size(original_tensor)
            @test size(cb_new.coeff_tensor) == expected_size
        end

        @testset "Error cases" begin
            atoms = [1, 2]
            ls = [1, 1]
            cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
            cb = cb_list[1]

            # Wrong length
            @test_throws ArgumentError reorder_atoms(cb, [1, 2, 3])
            @test_throws ArgumentError reorder_atoms(cb, [1])
        end

        @testset "Single atom case (N=1)" begin
            atoms = [1]
            ls = [2]
            cb_list = tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
            @test length(cb_list) > 0

            cb = cb_list[1]
            original_atoms = copy(cb.atoms)
            original_ls = copy(cb.ls)

            # For N=1, sorting doesn't change anything, but we can test the function
            new_atoms = [5]
            cb_new = reorder_atoms(cb, new_atoms)

            @test cb_new.atoms == [5]
            @test cb_new.ls == original_ls
            @test cb_new.Lf == cb.Lf
            @test cb_new.Lseq == cb.Lseq
            @test size(cb_new.coeff_tensor) == size(cb.coeff_tensor)
        end
    end

    # ------------------------------------------------------------------
    # Coverage uplift: surface area below was previously untested.
    # Tests assert each function's contract (documented show format,
    # isless total order, constructor error paths, sort-and-permute
    # invariant, field-equality round trip), not observed values.
    # ------------------------------------------------------------------

    @testset "Base.show(::CoupledBasis)" begin
        # Two-body, Lf = 0 — concrete instance via the canonical builder.
        cb = tesseral_coupled_bases_from_tesseral_bases([1, 1], [10, 20])[1]
        s = repr(cb)

        @test startswith(s, "CoupledBasis(")
        @test endswith(s, ")")
        # Every field documented in the source's show body appears.
        @test occursin("ls=", s)
        @test occursin("Lf=0", s)
        @test occursin("Lseq=", s)
        @test occursin("atoms=", s)
        # coeff_tensor is rendered as size(), not the full array.
        @test occursin("coeff_tensor=", s)
        @test occursin(string(size(cb.coeff_tensor)), s)
    end

    @testset "Base.isless: all-fields-equal returns false" begin
        # Two CoupledBasis values that agree on every field must satisfy
        # !isless(a, b) && !isless(b, a) (strict / irreflexive contract).
        # This hits the final `return false` path that is otherwise
        # unreachable when the values genuinely differ.
        cb1 = CoupledBasis([1, 1], 0, Int[], [10, 20], zeros(3, 3, 1))
        cb2 = CoupledBasis([1, 1], 0, Int[], [10, 20], zeros(3, 3, 1))
        @test !isless(cb1, cb2)
        @test !isless(cb2, cb1)
        @test !isless(cb1, cb1)
    end

    @testset "CoupledBasis_with_coefficient: full constructor" begin
        # N = 2 → Lseq length 0, tensor rank 3, Mf inferred from tensor's
        # last dim. Happy path first.
        ls = [1, 1]
        Lf = 0
        Lseq = Int[]
        atoms = [10, 20]
        coeff_tensor = zeros(3, 3, 1)
        coefficient = [0.5]
        multiplicity = 2
        clusters = [[10, 20]]

        cbc = CoupledBasis_with_coefficient(
            ls, Lf, Lseq, atoms, coeff_tensor, coefficient, multiplicity, clusters,
        )
        @test cbc.ls == ls
        @test cbc.Lf == Lf
        @test cbc.Lseq == Lseq
        @test cbc.atoms == atoms
        @test cbc.coeff_tensor == coeff_tensor
        @test cbc.coefficient == coefficient
        @test cbc.multiplicity == multiplicity
        @test cbc.clusters == clusters

        # Each ArgumentError site in the constructor body.
        # (1) Lseq length must be max(0, N-2). For N=2 that's 0.
        @test_throws ArgumentError CoupledBasis_with_coefficient(
            ls, Lf, [5], atoms, coeff_tensor, coefficient, multiplicity, clusters,
        )
        # (2) atoms length must equal length(ls).
        @test_throws ArgumentError CoupledBasis_with_coefficient(
            ls, Lf, Lseq, [10, 20, 30], coeff_tensor, coefficient, multiplicity, clusters,
        )
        # (3) ndims(coeff_tensor) must be length(ls) + 1.
        @test_throws ArgumentError CoupledBasis_with_coefficient(
            ls, Lf, Lseq, atoms, zeros(3, 3), coefficient, multiplicity, clusters,
        )
        # (4) length(coefficient) must equal size(coeff_tensor, R).
        @test_throws ArgumentError CoupledBasis_with_coefficient(
            ls, Lf, Lseq, atoms, coeff_tensor, [0.5, 0.5], multiplicity, clusters,
        )
        # (5) each clusters[i] must have length N.
        @test_throws ArgumentError CoupledBasis_with_coefficient(
            ls, Lf, Lseq, atoms, coeff_tensor, coefficient, multiplicity, [[10, 20, 30]],
        )
    end

    @testset "CoupledBasis_with_coefficient(cb, ...) convenience constructor" begin
        cb = tesseral_coupled_bases_from_tesseral_bases([1, 1], [10, 20])[1]
        Mf = size(cb.coeff_tensor, ndims(cb.coeff_tensor))
        coefficient = collect(Float64, 1:Mf)
        multiplicity = 7
        clusters = [cb.atoms]

        cbc = CoupledBasis_with_coefficient(cb, coefficient, multiplicity, clusters)

        # Five fields are forwarded from `cb`; the last three are the
        # explicit arguments.
        @test cbc.ls == cb.ls
        @test cbc.Lf == cb.Lf
        @test cbc.Lseq == cb.Lseq
        @test cbc.atoms == cb.atoms
        @test cbc.coeff_tensor == cb.coeff_tensor
        @test cbc.coefficient == coefficient
        @test cbc.multiplicity == multiplicity
        @test cbc.clusters == clusters
    end

    @testset "Base.show(::CoupledBasis_with_coefficient)" begin
        cb = tesseral_coupled_bases_from_tesseral_bases([1, 1], [10, 20])[1]
        Mf = size(cb.coeff_tensor, ndims(cb.coeff_tensor))
        cbc = CoupledBasis_with_coefficient(cb, fill(0.25, Mf), 3, [cb.atoms])
        s = repr(cbc)

        @test startswith(s, "CoupledBasis_with_coefficient(")
        @test endswith(s, ")")
        @test occursin("ls=", s)
        @test occursin("Lf=0", s)
        @test occursin("Lseq=", s)
        @test occursin("atoms=", s)
        @test occursin("coeff_tensor=", s)
        @test occursin(string(size(cbc.coeff_tensor)), s)
        @test occursin("coefficient=", s)
        @test occursin("multiplicity=3", s)
        @test occursin("clusters=", s)
    end

    @testset "reorder_atoms(::CoupledBasis_with_coefficient, ...)" begin
        # Mirror the existing two-body reorder_atoms test, additionally
        # pinning the coefficient-and-multiplicity invariant from the
        # function's docstring: those two fields ride the Mf dimension
        # (last axis) which is never permuted, so they pass through
        # unchanged. The clusters list, by contrast, lives in the per-site
        # axes and must be permuted by the same `p` that sorts the atoms.
        cb = tesseral_coupled_bases_from_tesseral_bases([3, 1], [1, 2])[1]
        Mf = size(cb.coeff_tensor, ndims(cb.coeff_tensor))
        coefficient = collect(Float64, 1:Mf)
        multiplicity = 11
        # Two-cluster orbit so the per-site permutation is observable.
        clusters = [[1, 2], [3, 4]]
        cbc = CoupledBasis_with_coefficient(cb, coefficient, multiplicity, clusters)
        orig_shape = size(cbc.coeff_tensor)

        cbc_new = reorder_atoms(cbc, [4, 1])

        # Atoms sorted, ls permuted to match.
        @test cbc_new.atoms == [1, 4]
        @test cbc_new.ls == [1, 3]
        @test cbc_new.Lf == cbc.Lf
        @test cbc_new.Lseq == cbc.Lseq
        # First N dims of coeff_tensor swap; the trailing Mf dim is fixed.
        @test size(cbc_new.coeff_tensor) ==
            (orig_shape[2], orig_shape[1], orig_shape[3])
        # Mf-dim invariants: coefficient and multiplicity must survive
        # the reorder unchanged.
        @test cbc_new.coefficient == coefficient
        @test cbc_new.multiplicity == multiplicity
        # Clusters live in the site axes: the permutation p = [2, 1]
        # (sorting [4, 1] → [1, 4]) flips the two columns of every cluster.
        @test cbc_new.clusters == [[2, 1], [4, 3]]

        # Error path: length(new_atoms) != length(ls).
        @test_throws ArgumentError reorder_atoms(cbc, [1, 2, 3])
        @test_throws ArgumentError reorder_atoms(cbc, [1])
    end

    @testset "convert_to_coupled_basis" begin
        cb = tesseral_coupled_bases_from_tesseral_bases([1, 1], [10, 20])[1]
        Mf = size(cb.coeff_tensor, ndims(cb.coeff_tensor))
        cbc = CoupledBasis_with_coefficient(cb, fill(0.5, Mf), 4, [cb.atoms])

        cb_back = convert_to_coupled_basis(cbc)

        # Result is a CoupledBasis and the five shared fields equal the
        # source.
        @test cb_back isa CoupledBasis
        @test cb_back.ls == cbc.ls
        @test cb_back.Lf == cbc.Lf
        @test cb_back.Lseq == cbc.Lseq
        @test cb_back.atoms == cbc.atoms
        @test cb_back.coeff_tensor == cbc.coeff_tensor

        # cb → cbc → cb_back round trip preserves the original five
        # fields (the coefficient/multiplicity are simply dropped).
        @test cb_back.ls == cb.ls
        @test cb_back.Lf == cb.Lf
        @test cb_back.Lseq == cb.Lseq
        @test cb_back.atoms == cb.atoms
        @test cb_back.coeff_tensor == cb.coeff_tensor
    end
end
