# `src/` refactoring review sweep

**Status**: complete (2026-06-01)

A multi-agent read-only review of all 23 `src/*.jl` files (~12.4k lines)
produced a graded findings list. The findings were worked one at a time
with per-item confirmation. Net result: a set of behavior-preserving
(bit-identical) cleanups landed across two tracks, with two items
deferred and two declined after investigation. This note is the
bird's-eye record; the architecture/performance subset has its own
canonical history under the spec folder.

## How the work was split

- **Standalone items on `main`** (bug fix, low-risk cleanup, dedup,
  API/format items): each a single confirmed commit.
- **Mid-sized architecture + hot-path subset** (categories A and C):
  isolated on a topic branch with a spec; merged via PR #33. The spec
  folder `docs/specs/260601-src-refactor/` is the canonical history for
  that subset (milestones, the A2 deferral reasoning, the bench numbers,
  and the Tier 2 four-axis panel outcome).

## Outcome by category

| Category | Scope | Outcome |
|---|---|---|
| Item 0 | `VaspIO` `parse_outcar` operator-precedence (`&&`/`\|\|`) | Fixed as a clarity change after verifying behavior-neutrality over 137,200 well-formed sequences (`f2db0ec`). Graded a "blocker" by the IO reviewer; the grade was wrong — the two parenthesizations diverge only on malformed non-numeric 7-column lines, where both error. |
| D | Docstrings / typos / named-tuple form; magic-number constants in `Clusters` (`DEFAULT_TOLERANCE`, `VIRTUAL_CELL_GRID_SIZE` / `NUM_VIRTUAL_CELLS`); `_`-prefix on ~15 internal helpers | Landed (`3020467`, `24eb843`, `5bf87b8`). |
| B3 | `SunnyExport` SALC-classification duplication | Extracted `_classify_salc` + `_unsupported_salc_string` (`8845e38`). |
| A1 | Angular-momentum coupling cache: hidden global + module-boundary breach | Added the `CoupledBases.cached_coupling_results` peek accessor; `SALCBases` no longer reaches into `_angular_momentum_cache` / the `:none` key (spec M1, `6a7e359`). |
| C1–C5 | SALC-build and torque SH-cache allocations | In-place projection write, buffer reuse in `is_translationally_equivalent_coupled_basis`, `(false, true)` tuple, and a combined SH value+gradient entry point `TesseralHarmonics.Zₗₘ_grad_unsafe` (one Legendre recursion instead of two; ~2.4× faster per-atom torque SH fill, bit-identical) (spec M2, `90f8755`). |

## E-items (API / format surface — each individually confirmed)

| Item | Decision |
|---|---|
| Drop `*_unsafe` from `TesseralHarmonics` `export` | Done (`a9855f6`). Explicit `using ..M: name` imports keep working for non-exported names, so only the one bare-`using` consumer (`Fitting`) needed an explicit import line; tests/benches were unaffected. |
| `calc_rmse` / `calc_r2score` → `_calc_rmse` / `_calc_r2score` | Done (`02514f5`); docstrings added; 8 qualified call sites in `Magesty.jl` updated. |
| `SALCBasis` `isotropy` positional → keyword | Done (`ae78baf`). The constructor's own docstring already documented it as a keyword, so this aligned code with documented API; no repo call site passed it positionally. |
| `kd_name` XML `<SpeciesOrder>` schema | **Declined.** The XML round trip is numerically exact (bit-identical, guarded by the existing save/load suite). On load, `kd_name` is rebuilt alphabetically (a documented, tested invariant) and `kd_int_list` is rebuilt consistently with it, so every consumer that resolves `kd_name[kd_int_list[atom]]` to an element string is order-independent. No species *index* is persisted or compared across the save/load boundary, so the lost user-defined ordering has no current consumer. A `<SpeciesOrder>` element would be an XML-format change touching the SCE-I/O linked site for no present numerical benefit. |

## Deferred / declined with reasoning

- **A2 — `Fitting` energy/torque kernel de-duplication: deferred.** The
  shared "skeleton" the review imagined does not cleanly exist. The
  energy element uses a deliberately different optimized structure
  (last-site split, scalar accumulation, `Z`-only) and cannot share the
  torque loop. The two torque accumulators are genuine near-twins and
  *could* be unified bit-identically, but one is on the package's hottest
  path (the threaded torque design-matrix build); routing it through a
  `view` adds indirection there for a maintainability-only gain (~30
  duplicated lines), which is not worth the regression risk. Revisit only
  if a future change already touches these kernels, or a benchmark shows
  the `view` path is neutral. (This was the performance reviewer's
  `[contention: performance]`.)
- **`kd_name` `<SpeciesOrder>`: declined** — see the E-items table.

## Lessons carried forward

- **Verify before accepting a reviewer's grade.** The "blocker" on the
  `VaspIO` finding and the "ordering bug" framing of `kd_name` both
  dissolved under direct verification; both turned out behavior-neutral.
  Reviewer severity is a prior, not a verdict.
- **Non-exported names are still importable.** Removing a name from
  `export` is not the breaking change an impact scan first suggested —
  explicit `using M: name` imports continue to resolve, so the blast
  radius is limited to bare-`using` consumers.
- **Mid-sized + multi-choice work earns a spec; the rest stays on
  `main`.** Splitting the architecture/hot-path subset onto a branch with
  a spec (PR #33) kept the standalone cleanups reviewable as small
  independent commits.
