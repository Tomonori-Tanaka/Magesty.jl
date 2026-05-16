# Design: `src/` layout refactor

## 1. 改名マッピング（確定版）

### 1.1 ファイル名統一（モジュール名と一致させる）

| 旧パス | 新パス | モジュール名 |
|---|---|---|
| `src/common/version.jl` | `src/Version.jl` | `Version` |
| `src/common/SortedCounter.jl` | `src/SortedCounters.jl` | `SortedCounters` |
| `src/utils/atomsbase_adapter.jl` | `src/AtomsBaseAdapter.jl` | `AtomsBaseAdapter` |
| `src/utils/xml_io.jl` | `src/XMLIO.jl` | `XMLIO` |

### 1.2 モジュール改称

| 旧名 | 新名 | 旧ファイル | 新ファイル |
|---|---|---|---|
| `MySphericalHarmonics` | **`TesseralHarmonics`** | `src/utils/MySphericalHarmonics.jl` | `src/TesseralHarmonics.jl` |
| `Basis` | **`CoupledBases`** | `src/types/Basis.jl` | `src/CoupledBases.jl` |
| `Optimize` | **`Fitting`** | `src/Optimize.jl` | `src/Fitting.jl` |

`TesseralHarmonics` は実球面調和（テッセラル形式）であることを学術用語で明示。
`CoupledBases` は SALCBases と命名が揃い、内容（CG 係数で組んだ coupled
basis）にも一致。`Fitting` は StatsAPI の `fit` 動詞とも対応する。

### 1.3 据え置き（変更なし）

主要モジュール: `Structures`, `Symmetries`, `Clusters`, `SALCBases`,
`SpinConfigs`, `AngularMomentumCoupling`, `RotationMatrix`, `ConfigParser`,
`SphericalHarmonicsTransforms`, `AtomCells`。

## 2. ディレクトリ構成（確定版: D-A フラット）

```
src/
  Magesty.jl                       # 主モジュール
  AngularMomentumCoupling.jl
  AtomCells.jl                     # 旧 types/
  AtomsBaseAdapter.jl              # 旧 utils/atomsbase_adapter.jl
  Clusters.jl
  ConfigParser.jl                  # 旧 utils/
  CoupledBases.jl                  # 旧 types/Basis.jl
  Fitting.jl                       # 旧 Optimize.jl
  RotationMatrix.jl                # 旧 utils/
  SALCBases.jl
  SortedCounters.jl                # 旧 common/SortedCounter.jl
  SphericalHarmonicsTransforms.jl  # 旧 utils/
  SpinConfigs.jl
  Structures.jl
  Symmetries.jl
  TesseralHarmonics.jl             # 旧 utils/MySphericalHarmonics.jl
  Version.jl                       # 旧 common/version.jl
  XMLIO.jl                         # 旧 utils/xml_io.jl
```

17 ファイル / 1 階層。`common/`, `types/`, `utils/` は廃止。

選定理由は plan ファイル参照（公開パッケージ調査の結果、同規模の
StatsAPI / Spglib / MultivariateStats / AtomsBase いずれもフラット～
1 階層浅い分割で運用されている）。

## 3. 影響範囲（参照箇所の調査結果）

### 3.1 `src/Magesty.jl`

- `include("common/version.jl")` → `include("Version.jl")`
- `include("common/SortedCounter.jl")` → `include("SortedCounters.jl")`
- `include("types/AtomCells.jl")` → `include("AtomCells.jl")`
- `include("types/Basis.jl")` → `include("CoupledBases.jl")`
- `include("utils/SphericalHarmonicsTransforms.jl")` →
  `include("SphericalHarmonicsTransforms.jl")`
- `include("utils/AngularMomentumCoupling.jl")` →
  `include("AngularMomentumCoupling.jl")`
- `include("utils/ConfigParser.jl")` → `include("ConfigParser.jl")`
- `include("utils/atomsbase_adapter.jl")` → `include("AtomsBaseAdapter.jl")`
- `include("utils/RotationMatrix.jl")` → `include("RotationMatrix.jl")`
- `include("utils/MySphericalHarmonics.jl")` →
  `include("TesseralHarmonics.jl")`
- `include("utils/xml_io.jl")` → `include("XMLIO.jl")`
- `include("Optimize.jl")` → `include("Fitting.jl")`
- 対応する `using .X` 群を `using .Basis` → `using .CoupledBases` 等に更新

### 3.2 src 内部の相互参照

- `Clusters.jl` L28: `using ..SortedCounters: SortedCounter` —
  パス変更なしで動く（モジュール名はそのまま）
- `SALCBases.jl` L15: `using ..SortedCounters: SortedCounter` — 同上
- `SALCBases.jl` 多数箇所: `Basis.CoupledBasis` → `CoupledBases.CoupledBasis`
- `SALCBases.jl`, `Fitting.jl`（旧 `Optimize.jl`）, ほか:
  `MySphericalHarmonics` 参照 → `TesseralHarmonics`
- `XMLIO.jl`: `using ..Version` — パス変更なし

### 3.3 test/

参照含む見込みファイル（grep で確認済み）:
- `test/runtests.jl`
- `test/component_test/test_MySphericalHarmonics.jl` →
  ファイル名も `test_TesseralHarmonics.jl` に
- `test/component_test/test_SortedCounter.jl` — モジュール名は既に
  `SortedCounters` なので中身次第（要 grep）
- `test/component_test/test_sphericart_agreement.jl`
- `test/examples/{chain,dimer}/test.jl`
- `bench/benchmark_sphericart.jl`, `bench/benchmark_spherical_harmonics.jl`

### 3.4 tools/

`src/` の改名・移動に追随する `include` / `using` の最小限の同期のみ。
構造変更や中身のリファクタは行わない（spec のスコープ外）。
着手時に `grep -rln "MySphericalHarmonics\|atomsbase_adapter\|xml_io\|
Optimize\|types/Basis" tools/` で確認。

### 3.5 docs / SPEC.md / CLAUDE.md

- `docs/src/api.md` の `@docs` ブロック:
  `MySphericalHarmonics.Zₗₘ` 等 → `TesseralHarmonics.Zₗₘ` 等
- `SPEC.md` のディレクトリ図と主要モジュール表
- `CLAUDE.md` 内の物理規約セクション「球面調和関数の規約」（旧名指し箇所）
- `docs/design-notes/` 内の旧名参照（履歴文書なので最小限のみ）

## 4. コミット粒度（実装順）

各コミット後に `make test-all` + `make test-tools` + `make build` が
緑であることを確認。各コミット独立に revert 可能。

### Commit 1 — ファイル名統一 + ディレクトリ平坦化

- ファイル移動 (`git mv`) 14 件:
  - `common/version.jl` → `Version.jl`
  - `common/SortedCounter.jl` → `SortedCounters.jl`
  - `types/AtomCells.jl` → `AtomCells.jl`
  - `types/Basis.jl` → `Basis.jl`（まだ改称しない）
  - `utils/SphericalHarmonicsTransforms.jl` → `SphericalHarmonicsTransforms.jl`
  - `utils/AngularMomentumCoupling.jl` → `AngularMomentumCoupling.jl`
  - `utils/ConfigParser.jl` → `ConfigParser.jl`
  - `utils/atomsbase_adapter.jl` → `AtomsBaseAdapter.jl`
  - `utils/RotationMatrix.jl` → `RotationMatrix.jl`
  - `utils/MySphericalHarmonics.jl` → `MySphericalHarmonics.jl`（まだ改称しない）
  - `utils/xml_io.jl` → `XMLIO.jl`
  - `Optimize.jl` (top-level, 変更なし)
- `src/Magesty.jl` の `include` パス更新（モジュール名は触らない）
- 空になった `common/`, `types/`, `utils/` ディレクトリ削除
- このコミット時点でモジュール名は据え置き → diff が最小・追跡しやすい

### Commit 2 — モジュール改称

- `MySphericalHarmonics.jl` → `TesseralHarmonics.jl` + `module MySphericalHarmonics` → `module TesseralHarmonics`
- `Basis.jl` → `CoupledBases.jl` + `module Basis` → `module CoupledBases`
- `Optimize.jl` → `Fitting.jl` + `module Optimize` → `module Fitting`
- 全 src 内 `using ..Basis` → `using ..CoupledBases` 等の参照書き換え
- `Basis.CoupledBasis` のドット参照 → `CoupledBases.CoupledBasis`
- test/ の参照書き換え + `test_MySphericalHarmonics.jl` →
  `test_TesseralHarmonics.jl`
- tools/ の参照書き換え（最小限）
- `docs/src/api.md` の `@docs` ブロック更新
- `Magesty.jl` 主モジュールの docstring 内の名指し更新
- `make test-all` + `make build` 緑確認

### Commit 3 — SPEC.md / CLAUDE.md / docs/design-notes 同期

- `SPEC.md` のディレクトリ図を新レイアウトに置き換え
- `CLAUDE.md` の物理規約セクション「球面調和関数の規約」内、
  `MySphericalHarmonics` 名指しを `TesseralHarmonics` に
- `docs/design-notes/` の名前参照を最小限で同期（履歴ノートなので
  あくまで「壊れた参照を直す」程度）
- `git grep -E "MySphericalHarmonics|module Basis\b|module Optimize\b|atomsbase_adapter|xml_io.jl|common/version|types/Basis|utils/MySphericalHarmonics"` でゼロを確認

## 5. リスクと緩和

- **見落とした参照**: 各コミット後に `git grep` で旧名・旧パスが
  ゼロであることを確認する。`make test-aqua` がモジュール構造の
  サニティチェックも兼ねる。
- **`@docs` の解決失敗**: Documenter は `@docs` ブロックでモジュール
  パスを解決する。Commit 2 で `docs/src/api.md` を更新し、`make build`
  が警告ゼロで通ることを確認。
- **`SortedCounters` の単複ずれ**: 既存の `module SortedCounters`（複数）
  をそのまま維持。中の型 `SortedCounter`（単数）も維持。Julia 標準
  ライブラリの `Strings` / `Tests` / `LinearAlgebra` と同じ規約。
- **`git mv` の履歴**: `git mv` で移動すれば `git log --follow` で
  履歴は追える。`git diff -M` も rename を検出する。
- **JET / Aqua の偽陽性**: モジュール改称後に未使用 export 等を Aqua
  が指摘する可能性。各コミット後に `make test-aqua` を実行。

## 6. 検証チェックリスト

各コミット後に以下を順に実行:

1. `make test-unit`
2. `make test-integration`
3. `make test-tools`
4. `make test-aqua`
5. `make test-jet`
6. `make build`（docs）
7. `git grep` で旧名・旧パスがゼロ:
   - `git grep -E "MySphericalHarmonics"` （Commit 2 後）
   - `git grep -E "module Basis\b"`（Commit 2 後）
   - `git grep -E "module Optimize\b"`（Commit 2 後）
   - `git grep -E "common/version|common/SortedCounter|types/AtomCells|types/Basis|utils/"` （Commit 1 後、Commit 2 でも）
