# 自作コンテナを DataStructures.jl で置き換える案

**Status**: **完了** (2026-05-14). CountingContainer は `OrderedDict` 置換、SortedContainer は新規 `SortedCounter` 型 + plain `Vector` + `sort!` once パターンに移行。

> **方針変更**: 当初本ノートは「`DataStructures.jl` の tree 系 (`SortedMultiSet` 等) に置き換える」案だったが、micro-benchmark の結果 `SortedMultiSet` が DataStructures.jl に存在しないこと、および tree 系は本ワークロードの規模では marginal 〜 逆効果であることが判明したため、spec (`docs/specs/260514-replace-sorted-container/`) で「plain container + sort once」パターンへ方針転換した。詳細は spec の design.md を参照。

**目的**: 自作の `SortedContainer.jl` (~360 行、3 コンテナ) と `CountingContainer.jl` (~120 行、1 コンテナ) を `DataStructures.jl` の成熟した実装に置き換え、保守コストとバグリスクを削減する。`DataStructures.jl` は既に依存。

## 動機

- 自作 4 型は `DataStructures.jl` にほぼ直接対応する型がある。
- `DataStructures.jl` は **既に Magesty.jl の依存** (`Basis.jl` の `counter` 関数で使用)。新規依存追加なし。
- 自作実装には現状もバグが残存している (refactor-sweep R8 Plan C 残課題: `SortedCountingUniqueVector::clear!` / `deleteat!` 未定義、`delete!` の O(N)、`copy` の無駄など)。
- Magesty.jl は pre-release で後方互換性の重荷がない → 置き換えの好機。

## 対応関係

### SortedContainer.jl

| Magesty 自作 | DataStructures.jl の対応物 |
|---|---|
| `SortedVector{T}` (重複許可) | `SortedMultiSet{T}` |
| `SortedUniqueVector{T}` | `SortedSet{T}` |
| `SortedCountingUniqueVector{T}` | `SortedMultiSet{T}` (multiplicity を内部追跡) |

`SortedMultiSet` の特徴:
- 平衡二分木 (AVL) ベース、挿入・削除・lookup すべて O(log N)
- `push!(s, x)` / `push!(s, x, n)` / `count(s, x)` / `delete!(s, x)` を提供
- iteration は sorted 順 (tree walk)

### CountingContainer.jl

**実施結果 (2026-05-14)**: `CountingUniqueVector{T}` → `OrderedDict{T, Int}` (DataStructures.jl) に置換し、`CountingContainer.jl` モジュールごと削除。

| Magesty 自作 | 採用した置き換え | 検討したが不採用 |
|---|---|---|
| `CountingUniqueVector{T}` (insertion order + count) | `OrderedDict{T, Int}` | `Accumulator{T, Int}` |

**`Accumulator` を採用しなかった理由**: 当初メモでは「`Clusters.jl` は順序に依存していない」と記載していたが、実装を再精査したところ誤り。`irreducible_clusters` (`src/Clusters.jl` L606-) は `cluster_dict[body][prim_atom_sc]` を iterate して、translationally equivalent なクラスタが既に追加済みなら `continue` する設計。「どのクラスタが canonical 代表として残るか」が iteration 順に依存しており、`Accumulator` (Dict ベース、順序保持なし) を使うと canonical 代表の atom 並びが変わって BasisSet 出力 (XML) が baseline と差分を起こす。

**`OrderedDict` で十分な根拠**: Clusters.jl で `CountingUniqueVector` から実利用しているのは 4 操作のみ:
- `CountingUniqueVector{Vector{Int}}()` → `OrderedDict{Vector{Int}, Int}()`
- `push!(cuv, val)` (重複時はカウント加算) → `od[k] = get(od, k, 0) + 1`
- `for x in cuv` (unique key を挿入順) → `for x in keys(od)`
- `getcounts(cuv, val)` → `cluster_dict[...][cluster]` (キー存在保証済み箇所では直接 lookup)

`AbstractVector` インタフェース (`length` / `size` / `getindex(::Int)` / `append!` / `==` / `isless` / `copy` 等) は Clusters.jl では一切呼ばれていなかった。`OrderedDict` は DataStructures.jl 提供で既に `using DataStructures` 済み、新規依存なし。

**連動箇所**: dead だった `using ..CountingContainer` を `BasisSets.jl` から削除。`Magesty.jl` の `include("common/CountingContainer.jl")`、`test/runtests.jl` の include と testset、`src/common/CountingContainer.jl`、`test/component_test/test_CountingContainer.jl` を全削除。

**数値結果**: 不変 (`test-unit` 6122/6122、`test-integration` 155/155、`test-jet`、`test-aqua` 全パス)。

## 性能 trade-off

| 観点 | 自作 (Vector ベース) | DataStructures.jl (Tree ベース) |
|---|---|---|
| 挿入 (sorted insertion) | O(N) (`insert!` で要素シフト) | O(log N) |
| 削除 | O(N) | O(log N) |
| Lookup | O(log N) (`searchsortedfirst`) | O(log N) |
| Iteration | O(N) cache-friendly (連続メモリ) | O(N) cache-unfriendly (tree walk) |
| メモリ | 連続 (少) | tree node 分散 (多) |
| 構築規模が小さい場合 | 有利 (定数オーバーヘッド小) | 不利 |
| 構築規模が大きい場合 | 不利 (O(N²)) | 有利 (O(N log N)) |

## 主要 caller

### SortedContainer.jl
- `BasisSets.jl`: `SortedCountingUniqueVector{Basis.CoupledBasis}` を coupled basis 構築で多用。
  - 構築時の push 回数 = 全 coupled basis の数 (N が大きいと O(N²) になりうる)。
  - 構築後の iteration が支配的なフェーズ (design matrix 構築) も多い。

### CountingContainer.jl
- `Clusters.jl` のみで使用 (`Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}` で interaction clusters を集計)。
  - 用途は単純な「unique 要素 + count」集計のみ。
  - 後段で `irreducible_cluster_dict::Dict{Int, SortedCountingUniqueVector{Vector{Int}}}` に変換 (こちらは sorted)。

## 着手前に必要なこと

### SortedContainer.jl
1. **ベンチマーク**: 既存の `BasisSet` 構築で Vector ベース vs `SortedMultiSet` の時間・メモリを測定。
   - 小規模 (dimer / chain): Vector 有利の可能性。
   - 大規模 (FeGe B20 2x2x2 / Fe BCC 2x2x2): Tree 有利の可能性。
   - 計測結果次第で「全置き換え」「型別置き換え」「現状維持」を判断。
2. **API gap の確認**: 自作型の `findall` / `deleteall!` 等が DataStructures.jl 側で同等に取れるか。
3. **callers の影響範囲**: `getcount` / `addcount!` 等の Magesty 固有 helper を全置換できるか。

### CountingContainer.jl
1. **caller 影響範囲**: `Clusters.jl` のみなので機械的に置換可能。
2. **ベンチマーク**: 内部実装はどちらも `Dict` ベースなので有意差は出にくいが、念のため `BasisSet` 構築時間で確認。

## 着手判断

### 優先度
- **CountingContainer.jl の置き換えが先**: caller が 1 ファイルに局在、性能 trade-off が小さい (どちらも `Dict` ベース)、順序保持の懸念なし。**圧倒的に簡単**。
- **SortedContainer.jl の置き換えは慎重に**: caller 多数、性能特性の trade-off (Vector vs Tree) があり要ベンチ。

### スコープ判断
- 単独 spec: 中規模リファクタ相当 → 着手する場合は `docs/specs/[YYMMDD]-replace-custom-containers/` を切って合意。
- 一体で進めるか分割するか: 上記の通り `CountingContainer` は単純、`SortedContainer` は要ベンチ。**分割推奨**。

### refactor-sweep との関係
- refactor-sweep R8 Plan C (`SortedContainer` のバグ修正) より優先度は低くない:
  - Plan C で bug を直してから置き換える: 安全だが二度手間。
  - 置き換えを先にやる: Plan C のバグは自動的に消える。
- 公開リリース前に決着させるのが望ましい (リリース後だと public API の変更を伴う恐れ)。

## やらないこと (現時点)

- 自作型を保持しつつ内部実装だけ DataStructures.jl に差し替える (薄いラッパー): 一見後方互換に見えるが、Magesty.jl はまだ public API として `SortedVector` 等を export していない (`SortedContainer` モジュール内の export のみ、再 export はしていない) ので、ラッパー化のメリットが小さい。完全置き換えで OK。

## 関連

- [refactor-sweep R8](refactor-sweep.md) — Plan B で重複集約は完了、Plan C で残バグ修正予定。本置き換え案を実施すると Plan C は不要になる。
- [SCE 公開 API 4 型構成](sce-public-api.md) — 公開リリース前の整理。本案も同期して進めるのが理想。
