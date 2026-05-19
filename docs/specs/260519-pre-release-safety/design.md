# Design: pre-release safety net

Status: draft (2026-05-19)

## Summary

Four small, independent changes share a common motivation — making
`SCEBasis` / `SCEFit` / `SCEModel` correctly self-describe what they
were built with, so that load-then-predict cannot silently use a
mismatched basis or a forgotten `torque_weight`, and so that the
unit-vector spin invariant is enforced rather than documented.

The basis fingerprint is the central design choice. The fingerprint is
a deterministic `UInt64` derived from the public-facing key-group
ordering of `salcbasis.salc_list` (the same ordering that the XML
serializer already preserves and that the design matrix relies on).
Anything that changes basis ordering — different `InteractionSpec`,
re-derived SALC for a different cutoff, accidentally swapped basis —
will produce a different fingerprint. `_check_basis` keeps its existing
`===` fast-path for in-session identity and adds a fingerprint
comparison for cross-process / cross-file cases. The fingerprint is
informational metadata only: predictions never branch on it.

`torque_weight` is added to `SCEModel` and to the model XML schema as a
plain field, mirroring `j0`. Older XML files default to `1.0`
(equivalent to the prior implicit behavior). This is additive and does
not break the XML round-trip invariant.

Spin unit-vector validation is a five-line change in the `SpinConfig`
inner constructor (per-column `isapprox(norm(c), 1.0; atol)`). The
tolerance is exposed as a keyword `atol_unit_norm` to keep test fixtures
that may have been written with slightly off-unit data working — but the
default is strict.

Missing tests are placed alongside the relevant existing tests, with no
new test infrastructure.

## Module layout

| Target | Change |
|---|---|
| `src/SALCBases.jl` | Add `salc_fingerprint(salcbasis::SALCBasis)::UInt64`; derive deterministically from the existing key-group ordering. |
| `src/Magesty.jl` (`SCEBasis` struct) | Add `salc_fingerprint::UInt64` as the 5th field. Provide an outer convenience constructor `SCEBasis(structure, symmetry, salcbasis, isotropy)` that computes the fingerprint via `salc_fingerprint(salcbasis)` and calls the default 5-arg constructor. The existing 3 call sites (`src/Magesty.jl:161, 1234, 1241`) keep their current 4-arg call shape and get the fingerprint for free. |
| `src/Magesty.jl` (`_check_basis`, lines 894-904) | The current body is a single `model.basis === dataset.basis \|\| throw(...)` expression. Restructure into two short-circuits: the existing identity check, then `model.basis.salc_fingerprint == dataset.basis.salc_fingerprint && return nothing`, with the throw at the end. Signature unchanged (`(model::SCEModel, dataset::SCEDataset)`); the `SCEFit` delegate at line 905 needs no change. |
| `src/Magesty.jl` (`SCEModel` struct, lines 265-278) | `SCEModel` currently has an **inner** constructor (lines 270-277) that validates `length(jphi) == length(basis.salcbasis.salc_list)`. Restructure: the inner constructor becomes the 4-arg form `SCEModel(basis, j0, jphi, torque_weight)` and retains the length check; add an outer 3-arg `SCEModel(basis, j0, jphi)` wrapper that defaults `torque_weight = 1.0` and delegates to the 4-arg inner. Existing 3-arg fixtures (e.g. `test/integration/square_lattice/test.jl:39`) keep working unchanged. |
| `src/SpinConfigs.jl` | Per-column unit-norm assertion in the inner constructor (`src/SpinConfigs.jl:101-157`); `atol_unit_norm` keyword (default `1e-6`, see Risks). |
| `src/XMLIO.jl` | Write / read `salc_fingerprint` attribute on the basis root and `torque_weight` element on the model root; tolerate absence on read (defaults: ignore disk value and recompute fingerprint for older basis XML; `torque_weight = 1.0` for older model XML). |
| `test/component/test_save_load.jl` | Round-trip test for both new fields; loadability of older XML synthesized on the fly with the new nodes removed. |
| `test/component/test_SpinConfigs.jl` | Reject a column with `norm = 1.01`; accept the same input under `atol_unit_norm = 1e-2`. |
| `test/component/test_fit_basis_check.jl` (new file) | Construct two `SCEModel` values whose `basis.salc_fingerprint` differs (via the default 5-arg constructor with a hand-picked fingerprint, since immutable structs cannot be mutated post-hoc), build an `SCEDataset` from one, and assert that `predict_energy(model_other, dataset)` raises the `_check_basis` error. |
| `test/integration/square_lattice/test.jl` | Add a new `@testset` that fits the existing isotropic basis (`fit(SCEFit, dataset, OLS())` against the analytic FM / AFM data) and verifies `predict_energy(SCEModel(fit), dataset)`, batched `predict_*` overloads, and XML round-trip of the fitted model. The existing direct-`SCEModel`-construction testset stays untouched. |
| pre-release cleanup design note | Remove the entry that owns the `SCEModel(fit)` round-trip; note in the "Issues found" section that the item moved into this spec. |
| `SPEC.md`, `docs/src/api.md`, `CHANGELOG.md` | Document new fields and the strict spin-norm invariant. |

## API

```julia
# SALCBases
function salc_fingerprint(b::SALCBasis)::UInt64

# SCEBasis: add field; outer constructor preserves the existing 4-arg shape
struct SCEBasis
    structure::Structure
    symmetry::Symmetry
    salcbasis::SALCBasis
    isotropy::Bool
    salc_fingerprint::UInt64   # NEW (5th field)
end

# Julia's implicit default constructor now takes all 5 positional fields:
#   SCEBasis(structure, symmetry, salcbasis, isotropy, salc_fingerprint)
# The new outer 4-arg constructor below computes the fingerprint and
# delegates to the 5-arg default form, so the existing 3 call sites at
# src/Magesty.jl:161, 1234, 1241 keep their 4-arg shape unchanged.
SCEBasis(structure::Structure, symmetry::Symmetry,
         salcbasis::SALCBasis, isotropy::Bool) =
    SCEBasis(structure, symmetry, salcbasis, isotropy,
             salc_fingerprint(salcbasis))

# SCEModel: add field. The existing inner constructor (which validates
# length(jphi) == n_salc) takes the new field; a 3-arg outer wrapper
# preserves the existing call shape.
struct SCEModel
    basis::SCEBasis
    j0::Float64
    jphi::Vector{Float64}
    torque_weight::Float64     # NEW; defaults to 1.0 in older XML and
                               # in the 3-arg convenience constructor

    function SCEModel(basis::SCEBasis, j0::Real,
                      jphi::AbstractVector{<:Real},
                      torque_weight::Real)
        n_salc = length(basis.salcbasis.salc_list)
        length(jphi) == n_salc || throw(ArgumentError(
            "length(jphi) ($(length(jphi))) must equal the number " *
            "of SALCs ($n_salc)"))
        return new(basis, j0, jphi, torque_weight)
    end
end

# Outer 3-arg wrapper: defaults torque_weight = 1.0
SCEModel(basis::SCEBasis, j0::Real, jphi::AbstractVector{<:Real}) =
    SCEModel(basis, j0, jphi, 1.0)

# SCEModel(f::SCEFit) carries torque_weight through
SCEModel(f::SCEFit) -> SCEModel  # propagates f.torque_weight

# SpinConfig inner ctor: new keyword
SpinConfig(energy, magmom_size, spin_directions, local_magfield;
           atol_unit_norm::Float64 = 1e-6)
```

`_check_basis` keeps its existing signature
(`(model::SCEModel, dataset::SCEDataset)` plus the `SCEFit` delegate)
and gains one extra short-circuit after the identity check:

```julia
function _check_basis(model::SCEModel, dataset::SCEDataset)
    model.basis === dataset.basis && return nothing
    model.basis.salc_fingerprint == dataset.basis.salc_fingerprint &&
        return nothing
    throw(ArgumentError(
        "the evaluation SCEDataset was built from a different SCEBasis " *
        "than the predictor ... " *
        "(salc_fingerprint $(model.basis.salc_fingerprint) vs " *
        "$(dataset.basis.salc_fingerprint))"))
end
_check_basis(f::SCEFit, dataset::SCEDataset) =
    _check_basis(SCEModel(f), dataset)  # unchanged
```

## Types and conventions

- `salc_fingerprint::UInt64` is in-session metadata only. It is not an
  input to any numerical computation. On XML load the on-disk value is
  ignored and the fingerprint is recomputed from the reconstructed
  `salcbasis`, so cross-Julia-version hash drift cannot cause spurious
  load-time mismatches. The on-disk value is written purely for
  external debugging / cross-save comparison.

- **Wire-tuple recipe** (authoritative). Iterate `salc_list ::
  Vector{Vector{CoupledBasis_with_coefficient}}` in the natural order
  (outer key-group index, then inner SALC index). Per element, fold in
  the structural identifiers of
  `CoupledBasis_with_coefficient` (`src/CoupledBases.jl:294-301`) that
  determine which design-matrix column the SALC produces:

  ```julia
  (coupled.ls,          # Vector{Int}: orbital ang-momenta per site
   coupled.Lf,          # Int: final coupled L
   coupled.Lseq,        # Vector{Int}: intermediate L sequence
   coupled.atoms,       # Vector{Int}: site indices
   coupled.multiplicity,
   length(coupled.coefficient))
  ```

  Numerical values of `coefficient` and `coeff_tensor` are deliberately
  **not** in the wire tuple: they are guaranteed by the XML round-trip
  to come back bit-exact, and the fingerprint's job is to guard SALC
  *ordering and identity*, not numerical drift. `multiplicity` is
  included because two coupled bases that agree on `(ls, Lf, Lseq,
  atoms)` but differ in `multiplicity` are distinct elements producing
  distinct design-matrix columns.

  The fingerprint is computed by feeding the per-element tuple through
  `Base.hash` and folding into a running `UInt64` accumulator. Concrete
  pseudo-code:

  ```julia
  function salc_fingerprint(b::SALCBasis)::UInt64
      h = hash(:Magesty_SALC_fingerprint_v1)
      for group in b.salc_list, coupled in group
          h = hash((coupled.ls, coupled.Lf, coupled.Lseq,
                    coupled.atoms, coupled.multiplicity,
                    length(coupled.coefficient)), h)
      end
      return h
  end
  ```

  The leading `:Magesty_SALC_fingerprint_v1` seed lets us bump the
  recipe (`_v2`) later without colliding with old saved values.
- `torque_weight` on `SCEModel` is the weighting that was active during
  the fit; it is **not** a knob for prediction. Prediction continues to
  ignore it. The field exists for reproducibility and audit only.
- Spin direction unit-norm constraint becomes a hard invariant at
  construction time. Existing in-memory fixtures that constructed
  `SpinConfig` directly are expected to comply; if any fail, they were
  bugs.

## Impact on linked sites

- [x] Spherical-harmonics convention (`TesseralHarmonics`): no change.
- [x] SCE coefficient XML (`save` / `load`): **additive**. New optional
      `salc_fingerprint` attribute on the basis root, new
      `torque_weight` element on the model root. Older XML without
      them must still load: basis loads recompute the fingerprint from
      the reconstructed SALC ordering; model loads default
      `torque_weight = 1.0`.
- [x] `Fitting` ↔ `SALCBasis`: fingerprint is derived from
      `SALCBasis.salc_list` key-group ordering; if that ordering ever
      changes deliberately, the fingerprint changes accordingly. No
      change to the design-matrix construction code.
- [ ] `.claude/agents/` references: nothing to update (no module renames).
- [x] `SPEC.md` / `docs/src/api.md`: add the new fields to the type
      diagrams; mention the new strict spin-norm invariant.

## Test strategy

- Component:
  - `test/component/test_save_load.jl`: round-trip of
    `salc_fingerprint` and `torque_weight`; loadability of an XML
    fixture lacking these new nodes (created on the fly by writing a
    file and then deleting the new nodes with `EzXML`).
  - `test/component/test_SpinConfigs.jl`: rejects a column with
    `norm = 1.01` under the default tolerance, and accepts the same
    input with `atol_unit_norm = 1e-2`.
  - `test/component/test_fit_basis_check.jl` (new file):
    construct a real `SCEBasis` `b1` via the convenience constructor,
    then build a second basis `b2` by calling the default 5-arg
    constructor with a forced fingerprint
    `forced_fp != salc_fingerprint(b1.salcbasis)`. Build an
    `SCEDataset` from `b1` and an `SCEModel` from `b2`; assert that
    `predict_energy(model_b2, dataset_b1)` throws an `ArgumentError`
    that mentions the fingerprint mismatch.
- Integration:
  - `test/integration/square_lattice/test.jl`: add a new `@testset`
    (alongside the existing isotropic testset) that calls
    `fit(SCEFit, dataset, OLS(); torque_weight = 0.5, verbosity = false)`
    on the analytic Heisenberg dataset, converts to `SCEModel(fit)`,
    and verifies (a) `predict_energy(SCEModel(fit), dataset)` agrees
    with `predict_energy(fit, dataset)` to machine precision,
    (b) the batched `predict_*` overloads over
    `AbstractVector{SpinConfig}` and `AbstractVector{<:AbstractMatrix}`
    agree with the dataset path, (c) `Magesty.save(SCEModel(fit), ...)`
    followed by `Magesty.load(SCEModel, ...)` preserves
    `torque_weight = 0.5` and `salc_fingerprint`. The existing
    direct-`SCEModel(basis, 0.0, jphi_list)` testset is left
    untouched.
- No new benchmarks needed; nothing in this spec touches a hot path.

## Risks and open items

- **Fingerprint stability**: The fingerprint must be reproducible
  *within a Julia version* so that two `SCEBasis` values constructed
  from the same SALC ordering in the same session compare equal.
  Cross-Julia-version stability is not required because we recompute
  the fingerprint on every XML load and ignore the on-disk value.
  `Base.hash` is sufficient under that constraint; we do not need a
  custom mixer. (The exact wire-tuple recipe is pinned in *Types and
  conventions* above.)

- **`SCEBasis` constructor cost on the XML load path**: adding the
  4-arg outer constructor that calls `salc_fingerprint(salcbasis)`
  means `Magesty.load(SCEBasis, path)` and
  `Magesty.load(SCEModel, path)` now do one extra pass over
  `salcbasis.salc_list`. The pass is `O(num_key_groups)` and touches
  only metadata (no SALC re-projection); empirically this is a
  fraction of the XML parse cost itself. We accept this; if profiling
  later shows it is non-trivial, we can add an internal 5-arg path
  that lets the loader skip the recompute when it trusts its own
  output.

- **Cross-basis mismatch test**: `SCEBasis` is immutable, so we
  cannot overwrite `salc_fingerprint` post-construction. The test
  builds the second basis by calling the default 5-arg constructor
  `SCEBasis(structure, symmetry, salcbasis, isotropy, forced_fp)`
  with a forced fingerprint that differs from
  `salc_fingerprint(salcbasis)`. This deliberately bypasses the
  convenience constructor to simulate "wrong fingerprint on disk"
  rather than building a second physical basis. The error path is
  thus exercised end-to-end through `predict_energy`. The default
  5-arg constructor is internal-only; we will not document or export
  it for external use.

- **Public exposure of `salc_fingerprint`**: the field is readable as
  `basis.salc_fingerprint`. We do not add a `salc_fingerprint(basis)`
  accessor function (only the `salc_fingerprint(::SALCBasis)`
  computation helper). Users who want to compare bases across two
  saves can read the field directly.

- **`atol_unit_norm` default**: set to `1e-6` rather than `1e-8`.
  Spin direction columns in EMBSET and TOML inputs may have been
  rounded to ~6-7 significant figures by an upstream writer, and
  `1e-8` could falsely reject legitimate inputs. `1e-6` still catches
  any meaningful violation (a 1e-4 deviation already shifts `Zlm`
  noticeably). Tighten later if real fixtures show no false negatives.

- **`torque_weight` semantics when an `SCEModel` is built outside
  `fit`**: the 3-arg `SCEModel(basis, j0, jphi)` convenience
  constructor defaults `torque_weight = 1.0`, matching prior implicit
  behavior. This is what `test/integration/square_lattice/test.jl:39`
  uses today and we deliberately keep it working.
