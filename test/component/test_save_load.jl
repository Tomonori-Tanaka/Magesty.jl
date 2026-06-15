using Test
using Magesty

# Fixtures -------------------------------------------------------------

const SL_DIMER_TOML = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")
const SL_FEPT_TOML =
    joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2", "input.toml")
const SL_FEPT_EMBSET =
    joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2", "EMBSET")
const SL_FEGE_TOML =
    joinpath(@__DIR__, "..", "integration", "fege_2x2x2", "input.toml")
const SL_BASELINE_DIR = joinpath(@__DIR__, "baselines")


# Structural + numerical equality of two SALC bases. After an XML
# round-trip the SALC basis is a fresh object (`===` does not hold), so
# the design-matrix `(l, m, site)` ordering is verified by comparing the
# salc_list contents element-by-element. The default `atol = 1e-10`
# matches the SALC writer precision; relax to `1e-8` for cross-platform
# fresh-build comparisons where BLAS / LAPACK round-off accumulates.
function _assert_salcbasis_equal(a, b; atol::Real = 1e-10)
    @test length(a.salc_list) == length(b.salc_list)
    for (group_a, group_b) in zip(a.salc_list, b.salc_list)
        @test length(group_a) == length(group_b)
        for (cbc_a, cbc_b) in zip(group_a, group_b)
            @test cbc_a.ls == cbc_b.ls
            @test cbc_a.Lf == cbc_b.Lf
            @test cbc_a.Lseq == cbc_b.Lseq
            @test cbc_a.atoms == cbc_b.atoms
            @test cbc_a.multiplicity == cbc_b.multiplicity
            @test isapprox(cbc_a.coefficient, cbc_b.coefficient, atol = atol)
            @test size(cbc_a.coeff_tensor) == size(cbc_b.coeff_tensor)
            @test isapprox(cbc_a.coeff_tensor, cbc_b.coeff_tensor, atol = atol)
        end
    end
    return nothing
end

# The XML writer stamps the running package version into the <version>
# element. That field is environment metadata, not part of the schema or
# the SALC coefficients this regression guards, so mask it before the
# byte-equality comparison to keep the test stable across version bumps.
_mask_version(xml::AbstractString) =
    replace(xml, r"<version>.*?</version>" => "<version></version>")

# Tests ----------------------------------------------------------------

@testset "save / load" begin
    @testset "extension dispatch" begin
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        @test_throws ArgumentError Magesty.save(basis, "model.txt")
        @test_throws ArgumentError Magesty.save(basis, "model")
        @test_throws ArgumentError Magesty.load(SCEBasis, "basis.json")
        @test_throws ArgumentError Magesty.load(SCEModel, "model.yaml")

        # An uppercase `.XML` extension is accepted.
        mktempdir() do dir
            path = joinpath(dir, "BASIS.XML")
            Magesty.save(basis, path)
            @test isfile(path)
        end
    end

    @testset "SCEBasis round-trip" begin
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        mktempdir() do dir
            path = joinpath(dir, "basis.xml")
            Magesty.save(basis, path)
            reloaded = Magesty.load(SCEBasis, path)

            @test reloaded isa SCEBasis
            @test reloaded.isotropy == basis.isotropy
            @test reloaded.symmetry.tol == basis.symmetry.tol
            @test reloaded.structure.supercell.num_atoms ==
                  basis.structure.supercell.num_atoms
            @test reloaded.structure.supercell.lattice_vectors ≈
                  basis.structure.supercell.lattice_vectors
            @test reloaded.structure.supercell.x_frac ≈ basis.structure.supercell.x_frac
            _assert_salcbasis_equal(basis.salcbasis, reloaded.salcbasis)
        end
    end

    @testset "SCEBasis round-trip preserves asymmetric lattice" begin
        # Regression: the lattice writer emitted matrix rows under the
        # `a$i` tags while the reader parses each `a$i` into column `i`,
        # so a save / load round trip transposed `lattice_vectors`. A
        # symmetric lattice (`L == Lᵀ`) hides the bug; an asymmetric
        # (here upper-triangular) cell exposes it. Reciprocal vectors,
        # Cartesian rotations, the SALC projection, and the design matrix
        # are all rebuilt from the reloaded lattice, so a transpose
        # silently corrupts every low-symmetry cell.
        basis = SCEBasis(;
            lattice = [3.0 1.0 0.5; 0.0 3.0 0.7; 0.0 0.0 3.0],
            kd = [:X], kd_list = [1, 1],
            positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1]],
            periodicity = (false, false, false),
            interaction = (body1 = (lmax = Dict(:X => 0),),
                           body2 = (lsum = 2, cutoff = Dict((:X, :X) => -1.0))),
            isotropy = true, verbosity = false,
        )
        L = basis.structure.supercell.lattice_vectors
        # The fixture must be asymmetric, or the round trip could not
        # distinguish `L` from `Lᵀ`.
        @test L != permutedims(L)
        mktempdir() do dir
            path = joinpath(dir, "basis.xml")
            Magesty.save(basis, path)
            reloaded = Magesty.load(SCEBasis, path)
            L_reloaded = reloaded.structure.supercell.lattice_vectors
            @test L_reloaded ≈ L
            @test !isapprox(L_reloaded, permutedims(L))
        end
    end

    @testset "SCEModel round-trip" begin
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        n_salc = length(basis.salcbasis.salc_list)
        # Construct a model directly from synthetic coefficients — the
        # round-trip is about serialization, not about the fit.
        model = SCEModel(basis, -3.141592653589793, collect(1:n_salc) ./ 7)
        mktempdir() do dir
            path = joinpath(dir, "model.xml")
            Magesty.save(model, path)
            reloaded = Magesty.load(SCEModel, path)

            @test reloaded isa SCEModel
            @test reloaded.j0 ≈ model.j0
            @test reloaded.jphi ≈ model.jphi
            @test reloaded.basis.isotropy == basis.isotropy
            @test reloaded.basis.symmetry.tol == basis.symmetry.tol
            _assert_salcbasis_equal(model.basis.salcbasis, reloaded.basis.salcbasis)
        end
    end

    @testset "SCEModel j0 / jphi exact round-trip" begin
        # The XML writer must preserve every Float64 mantissa bit in
        # `j0` and `jphi`, so that save -> load reproduces the input
        # values bit-for-bit (Float64 round-trip). The previous
        # `%.15e` format lost 1 ulp on adversarial inputs (e.g.
        # `nextfloat(1.0)` formatted as `1.000000000000000e+00`).
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        n_salc = length(basis.salcbasis.salc_list)
        # Pick values that exercise the worst case for fixed-precision
        # decimal: ulp-adjacent floats and irrationals are the ones
        # that `%.15e` used to mangle. The vector mixes them across
        # every SALC index so no single Magesty-side rounding shortcut
        # can mask a regression.
        adversarial_floats = Float64[
            nextfloat(1.0),                  # 1 ulp above 1
            prevfloat(1.0),                  # 1 ulp below 1
            nextfloat(0.1),
            1.2345678901234567e-15,
            -3.7e-200,
            1.0 + eps(Float64),
            Float64(pi),
            sqrt(2.0),
        ]
        jphi = Float64[
            adversarial_floats[mod1(i, length(adversarial_floats))]
            for i in 1:n_salc
        ]
        j0 = nextfloat(-3.14159265358979)
        model = SCEModel(basis, j0, jphi)
        mktempdir() do dir
            path = joinpath(dir, "model.xml")
            Magesty.save(model, path)
            reloaded = Magesty.load(SCEModel, path)

            @test reloaded.j0 === model.j0       # IEEE-754 bit equality
            @test reloaded.jphi == model.jphi    # element-wise equality
        end
    end

    @testset "save(::SCEFit) delegates to SCEModel" begin
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        SC = Magesty.SpinConfigs.SpinConfig
        configs = [
            SC(-1.0, [1.0, 1.0],
               [0.0 0.0; 0.0 0.0; 1.0 1.0],
               [0.1 0.0; 0.0 -0.1; 0.0 0.0]),
            SC(1.0, [1.0, 1.0],
               [0.0 0.0; 0.0 0.0; -1.0 1.0],
               [-0.1 0.0; 0.0 0.1; 0.0 0.0]),
        ]
        dataset = SCEDataset(basis, configs)
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3, verbosity = false)

        # save(::SCEFit) is a convenience wrapper around save(SCEModel(f)):
        # the two files must be byte-identical.
        mktempdir() do dir
            path_fit = joinpath(dir, "from_fit.xml")
            path_model = joinpath(dir, "from_model.xml")
            Magesty.save(f, path_fit)
            Magesty.save(SCEModel(f), path_model)
            @test read(path_fit) == read(path_model)

            reloaded = Magesty.load(SCEModel, path_fit)
            @test reloaded.j0 === f.j0
            @test reloaded.jphi == f.jphi
        end
    end

    @testset "load SCEBasis from an SCEModel XML" begin
        basis = SCEBasis(SL_DIMER_TOML; verbosity = false)
        n_salc = length(basis.salcbasis.salc_list)
        model = SCEModel(basis, 1.5, collect(1.0:n_salc))
        mktempdir() do dir
            model_path = joinpath(dir, "model.xml")
            Magesty.save(model, model_path)

            # Magesty.load(SCEBasis, ...) accepts an SCEModel XML, dropping <JPhi>.
            from_model = Magesty.load(SCEBasis, model_path)
            @test from_model isa SCEBasis
            _assert_salcbasis_equal(basis.salcbasis, from_model.salcbasis)

            # Magesty.load(SCEModel, ...) on a basis-only XML is a clear error.
            basis_path = joinpath(dir, "basis.xml")
            Magesty.save(basis, basis_path)
            @test_throws ArgumentError Magesty.load(SCEModel, basis_path)
        end
    end

    @testset "basis identity (design matrices) after reload" begin
        # The real round-trip guarantee: a dataset built from the reloaded
        # basis produces the same design matrices as the original, i.e.
        # the `(l, m, site)` ordering survives serialization. Since the
        # SALC coupling tensors are now serialized at `%.17e` (the
        # Float64-round-trip-safe precision), equality is bit-exact, not
        # just `≈`.
        basis = SCEBasis(SL_FEPT_TOML; verbosity = false)
        dataset = SCEDataset(basis, SL_FEPT_EMBSET)
        mktempdir() do dir
            path = joinpath(dir, "fept_basis.xml")
            Magesty.save(basis, path)
            reloaded = Magesty.load(SCEBasis, path)
            dataset_reloaded = SCEDataset(reloaded, SL_FEPT_EMBSET)

            @test dataset_reloaded.X_E == dataset.X_E
            @test dataset_reloaded.X_T == dataset.X_T
            # Observations come straight from the EMBSET file — exact.
            @test dataset_reloaded.y_E == dataset.y_E
            @test dataset_reloaded.y_T == dataset.y_T
        end
    end

    # Regression against committed baselines: catches accidental schema
    # drift in the XML writer / reader and unintended changes to SALC
    # coefficients. The load -> re-save loop is fully deterministic
    # (byte equality). The fresh-build loop additionally exercises the
    # eigensolver-based SALC construction, where the canonical gauge fix
    # in `_compute_salc_groups` collapses the degenerate-eigenspace
    # ambiguity to numerical noise — but BLAS / LAPACK round-off
    # (different SIMD lane widths on x86 vs arm64) still accumulates
    # below ~1e-8, so the fresh-build check uses `≈` at that tolerance
    # instead of byte equality.
    @testset "regression vs committed baseline" begin
        @testset "load -> re-save (byte equality)" begin
            for (name, T) in (
                ("fept_basis", SCEBasis),
                ("fept_model", SCEModel),
                ("fege_basis", SCEBasis),
            )
                baseline = joinpath(SL_BASELINE_DIR, "$name.xml")
                obj = Magesty.load(T, baseline)
                mktempdir() do dir
                    path = joinpath(dir, "$name.xml")
                    Magesty.save(obj, path)
                    @test _mask_version(read(path, String)) ==
                          _mask_version(read(baseline, String))
                end
            end
        end

        @testset "fresh build ≈ committed baseline (SCEBasis, atol=1e-8)" begin
            for (name, toml) in (
                ("fept_basis", SL_FEPT_TOML),
                ("fege_basis", SL_FEGE_TOML),
            )
                baseline = Magesty.load(SCEBasis, joinpath(SL_BASELINE_DIR, "$name.xml"))
                fresh = SCEBasis(toml; verbosity = false)
                _assert_salcbasis_equal(baseline.salcbasis, fresh.salcbasis; atol = 1e-8)
            end
        end
    end
end
