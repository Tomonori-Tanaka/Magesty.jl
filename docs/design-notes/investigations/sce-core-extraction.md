# Extracting a shared SCE-model core package

**Status**: complete — investigation only, no action taken (2026-06-14)

Conclusion: a shared "represent + serialize + evaluate a finished SCE
model" tier is cleanly separable from the "construct / fit" tier, and is
worth extracting eventually. Not done now: both packages are at `0.1.0`
and the on-disk format / SALC ordering are still moving, so the lockstep
coordination cost of a split outweighs the benefit today. Recorded as an
on-hold idea. When revisited, prefer an in-repo subpackage over a third
standalone repository, and fix the Spglib-on-load wart first.

## Background

SCE models built by Magesty are consumed downstream by the Monte Carlo
package `SpinClusterMC.jl`, with active-learning and Bayesian-modeling
consumers anticipated. The question: should the SCE model struct + XML
I/O (and the minimal evaluation surface) be factored into a shared
package that all consumers depend on, rather than each depending on the
full Magesty stack?

## What the dependency map shows

The construct/fit boundary and the represent/serialize/evaluate boundary
already coincide in practice.

- **Downstream surface is tiny.** `SpinClusterMC` touches exactly three
  Magesty symbols: the `CoupledBasis_with_coefficient` data type, the
  `read_salcbasis_from_xml` parse, and `build_all_real_bases` (tesseral
  CG enumeration, used only by its reference evaluation path to build a
  CG table). It does not use fitting, the design-matrix builders,
  optimization, DFT/input parsing, or Magesty's tesseral-harmonics
  routines (it uses SpheriCart directly).
- **Coupling is load-time only.** The MC sweep never calls back into
  Magesty; all Magesty data is unpacked into local structs at load.
  No version-pin hacks, no vendored copies, no type piracy.
- **Candidate core tier**: the `CoupledBases` / `SALCBases` data types,
  `Structure` / `Symmetry`, XML I/O, tesseral harmonics, angular-momentum
  coupling, and the `predict_energy` / `predict_torque` evaluation.
  Dependencies: `DataStructures`, `StaticArrays`, `LegendrePolynomials`,
  `EzXML`, `WignerSymbols`, `Combinat`.
- **Stays in Magesty (construct/fit)**: the design-matrix builders and
  regression, cluster generation, the SALC-basis constructor, and the
  input-spec / AtomsBase adapters. These pull the heavy deps
  (`GLMNet`, `Roots`, `AtomsBase`, `Unitful`) that downstream consumers
  do not want transitively.

## Frictions that block a clean leaf package

1. **Spglib-on-load.** Loading recomputes the `Symmetry` object via
   Spglib because the XML stores only `tolerance_sym`, not the full
   symmetry data. So the load path is not pure data and would drag
   Spglib into the core. Serializing the complete symmetry into the XML
   removes this — a format change (linked site: XML I/O round trip and
   SALC ordering) that needs the round-trip agreement test, and is
   worthwhile on its own merits.
2. **Duplicated XML parsing.** `SpinClusterMC` re-implements XML parsing
   in its reference path for cross-checking. That independent parser is
   a verification asset, but it also means the on-disk format already
   has a second reader that can drift from Magesty's.

## Why not act now

Both packages are `0.1.0` and the format / SALC ordering are the
volatile-yet-critical linked sites. A separate repository forces a
three-way lockstep version bump on every format change during exactly
the phase when the format is most likely to move.

## If revisited

- Prefer an **in-repo subpackage** (mirroring the prior CLI extraction
  into its own package directory) over a third standalone repo: Magesty
  tracks it by path dependency, downstream depends on the thin core
  without the heavy fitting deps, and there is no independent release
  cadence to coordinate until one is genuinely needed.
- Do the de-risking work first regardless of the split decision:
  (a) serialize full symmetry into the XML to make the load path pure
  data and drop Spglib from the core; (b) firewall the core modules so
  they have no upward dependency on the fitting / construction code, and
  enforce it with a dependency check.
- Treat the XML round-trip / agreement test as the shared contract; it
  becomes the core package's test suite.
- Natural trigger to start: a second real consumer (active learning)
  materializing, or the on-disk format stabilizing.
