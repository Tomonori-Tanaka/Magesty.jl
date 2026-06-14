using TOML
using Printf
using Test
using LinearAlgebra
using Magesty


@testset "Square Lattice Tests" begin

    spin_directions_fm = [
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0;
        1.0 1.0 1.0 1.0]
    spin_directions_afm = [
        0.0  0.0  0.0 0.0;
        0.0  0.0  0.0 0.0;
        1.0 -1.0 -1.0 1.0]

    @testset "Hypothetical SCE model (isotropic)" begin
        input = TOML.parse(open(joinpath(@__DIR__, "input_isotropic.toml"), "r"))
        basis = SCEBasis(input; verbosity = false)
        Magesty.save(basis, joinpath(@__DIR__, "system_isotropic.xml"))
        # Heisenberg model: E = Σ_{<ij>} Jij * Si·Sj, Jij = -1.0 eV (ferromagnetic) for all pairs.
        #
        # Square lattice (4 atoms, PBC): 8 NN pairs + 8 diagonal (NNN) pairs within the supercell.
        # Each SALC has a symmetry coefficient c = norm of the first basis's coefficient vector.
        # For Lf=0, the coefficient vector has length 1, so c = coefficient[1].
        # For Lf=2, the vector has length 5 (one per Mf), so we take the norm.
        #
        # Conversion from physical Jij to SCE coefficient jphi:
        #   jphi = Jij / (c × √3)
        # where the √3 comes from the l=1 spherical harmonics coupling:
        #   (4π) × (1/√3) × (3/4π) = √3  per pair.
        J = -1.0
        # Read SALC coefficients from the basis to be robust against ordering changes.
        # When a SALC has multiple bases, the norm of the first basis's coefficient vector is used.
        salc_coeffs = [norm(salc[1].coefficient) for salc in basis.salcbasis.salc_list]
        jphi_list = [J / (c * sqrt(3)) for c in salc_coeffs]
        model = SCEModel(basis, 0.0, jphi_list)

        # Ferromagnetic configuration: all spins along +z, Si·Sj = 1 for all pairs.
        # E = Jij × (8 NN + 8 diagonal) × 1 = -1.0 × 16 = -16 eV
        energy_analytic_fm = J * 8 + J * 8  # -16.0 eV
        @test predict_energy(model, spin_directions_fm) ≈ energy_analytic_fm atol = 1e-6

        # Antiferromagnetic configuration: S1=S4=+z, S2=S3=-z.
        # NN pairs:      all Si·Sj = -1  → Σ_{NN} = -8
        # Diagonal pairs: (1,4) and (2,3) are same spin → Si·Sj = +1 → Σ_{diag} = +8
        # E = J × (-8 + 8) = 0 eV
        energy_analytic_afm = 0.0
        @test predict_energy(model, spin_directions_afm) ≈ energy_analytic_afm atol = 1e-6
    end

    @testset "SCEFit / SCEModel round trip" begin
        # Exercises the production path SCEDataset → SCEFit → SCEModel
        # (vs. the direct-SCEModel construction in the testset above)
        # together with the batched / single-config predict_* overloads
        # and the save/load round trip routed through the recomputed
        # basis fingerprint.
        input = TOML.parse(open(joinpath(@__DIR__, "input_isotropic.toml"), "r"))
        basis = SCEBasis(input; verbosity = false)

        # Heisenberg energy / local_magfield on the 1×1 PBC supercell.
        # Within the cell, NN bonds {(1,2),(1,3),(2,4),(3,4)} each appear
        # twice across PBC and diagonal bonds {(1,4),(2,3)} each appear
        # four times — matching the 8 NN + 8 diagonal counts the
        # isotropic testset uses. local_magfield[i] = -m_i · Σ_j J_ij S_j
        # with m_i = 1; the multiplicity is folded into the per-neighbor
        # weight.
        J = -1.0
        nn_pairs   = ((1, 2), (1, 3), (2, 4), (3, 4))
        diag_pairs = ((1, 4), (2, 3))
        nbr_table  = (
            ((2, 2.0), (3, 2.0), (4, 4.0)),
            ((1, 2.0), (3, 4.0), (4, 2.0)),
            ((1, 2.0), (2, 4.0), (4, 2.0)),
            ((1, 4.0), (2, 2.0), (3, 2.0)),
        )

        function make_config(S::AbstractMatrix{Float64})
            E = 0.0
            for (i, j) in nn_pairs
                E += 2.0 * J * dot(@view(S[:, i]), @view(S[:, j]))
            end
            for (i, j) in diag_pairs
                E += 4.0 * J * dot(@view(S[:, i]), @view(S[:, j]))
            end
            H = zeros(3, 4)
            for i in 1:4
                for (j, mult) in nbr_table[i]
                    H[:, i] .+= (-1.0 * J * mult) .* @view(S[:, j])
                end
            end
            return Magesty.SpinConfigs.SpinConfig(E, ones(4), copy(S), H)
        end

        configs = Magesty.SpinConfigs.SpinConfig[
            # FM: all +z.
            make_config([0.0  0.0  0.0  0.0;
                         0.0  0.0  0.0  0.0;
                         1.0  1.0  1.0  1.0]),
            # AFM checkerboard: (1,4)=+z, (2,3)=-z.
            make_config([0.0  0.0  0.0  0.0;
                         0.0  0.0  0.0  0.0;
                         1.0 -1.0 -1.0  1.0]),
            # Single flip from FM: atom 1 down, rest up.
            make_config([ 0.0  0.0  0.0  0.0;
                          0.0  0.0  0.0  0.0;
                         -1.0  1.0  1.0  1.0]),
            # Row stripe: atoms (1,2) up, (3,4) down.
            make_config([0.0  0.0  0.0  0.0;
                         0.0  0.0  0.0  0.0;
                         1.0  1.0 -1.0 -1.0]),
            # Mixed: atoms (1,2) along +x, atoms (3,4) along +z.
            make_config([1.0  1.0  0.0  0.0;
                         0.0  0.0  0.0  0.0;
                         0.0  0.0  1.0  1.0]),
            # Triaxial: each atom on a different Cartesian axis pair.
            make_config([1.0  0.0  0.0  0.0;
                         0.0  0.0  1.0  0.0;
                         0.0 -1.0  0.0  1.0]),
            # Twisted: four distinct unit directions, all non-trivial.
            make_config(hcat(
                normalize([1.0, 0.0, 0.0]),
                normalize([0.0, 1.0, 0.0]),
                normalize([1.0, 1.0, 0.0]),
                normalize([1.0, 0.0, 1.0]),
            )),
        ]

        dataset = SCEDataset(basis, configs)
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.5, verbosity = false)
        model = SCEModel(f)

        # (c) SCEFit and SCEModel produce identical predictions.
        E_dataset = predict_energy(model, dataset)
        T_dataset = predict_torque(model, dataset)
        @test E_dataset ≈ predict_energy(f, dataset)
        T_from_fit = predict_torque(f, dataset)
        @test all(T_dataset[i] ≈ T_from_fit[i] for i in eachindex(T_dataset))

        # (d) Batched predict_* overloads over Vector{SpinConfig} and
        # Vector{<:AbstractMatrix} agree with the dataset path.
        sd_list = [c.spin_directions for c in configs]
        @test predict_energy(model, configs) ≈ E_dataset
        @test predict_energy(model, sd_list) ≈ E_dataset
        T_configs = predict_torque(model, configs)
        T_sd = predict_torque(model, sd_list)
        @test all(T_configs[i] ≈ T_dataset[i] for i in eachindex(T_dataset))
        @test all(T_sd[i] ≈ T_dataset[i] for i in eachindex(T_dataset))

        # (d') Single-config predict_* paths agree with the dataset path
        # entry by entry, both via SCEModel and via SCEFit, and for both
        # SpinConfig and bare spin-direction matrices. The explicit
        # atol covers configs whose predicted energy or torque rounds to
        # ~0 (e.g. the AFM checkerboard): there the dataset-path
        # evaluation (X_E * jphi) and the per-config basis evaluation
        # arrive at the same value through different summation orders,
        # and cross-BLAS rounding (Accelerate vs OpenBLAS) leaves a
        # residual on the order of eps(Float64); the default rtol-only
        # ≈ degenerates at zero. 1e-10 eV is six orders of magnitude
        # below the SCE energy scale and well above the ULP noise.
        for (i, sc) in enumerate(configs)
            @test predict_energy(model, sc) ≈ E_dataset[i] atol = 1e-10
            @test predict_energy(model, sc.spin_directions) ≈ E_dataset[i] atol = 1e-10
            @test predict_energy(f, sc) ≈ E_dataset[i] atol = 1e-10
            @test predict_torque(model, sc) ≈ T_dataset[i] atol = 1e-10
            @test predict_torque(model, sc.spin_directions) ≈ T_dataset[i] atol = 1e-10
            @test predict_torque(f, sc) ≈ T_dataset[i] atol = 1e-10
        end

        # (e) save / load round trip: the reloaded basis recomputes the
        # same fingerprint, so _check_basis accepts the in-memory dataset
        # and predictions reproduce to machine precision.
        mktempdir() do dir
            path = joinpath(dir, "square_model.xml")
            Magesty.save(model, path)
            reloaded = Magesty.load(SCEModel, path)
            @test reloaded.basis.salc_fingerprint == model.basis.salc_fingerprint
            @test predict_energy(reloaded, dataset) ≈ E_dataset
            T_reloaded = predict_torque(reloaded, dataset)
            @test all(T_reloaded[i] ≈ T_dataset[i] for i in eachindex(T_dataset))
        end
    end

end
