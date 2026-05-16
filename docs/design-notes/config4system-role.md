# Config4System の役割の見直し

**Status**: spec 化に進む（2026-05-14 提案 / 2026-05-16 案 B 深掘り）。
spec: [`docs/specs/260516-typed-input-spec/`](../specs/260516-typed-input-spec/)。
トリガーになったのはワイルドカード species 指定（`*-*`, `Fe-*`, `Fe-Fe`）の
新規導入。dict パース層に substantive な責務が増えるため、案 A での後付けより
案 B 移行と一度でやるほうが筋がよい。

**関連**: [`docs/specs/260514-sce-public-api/`](../specs/260514-sce-public-api/)
（SCE 公開 API リファクタ。merged）

## 背景

SCE 公開 API リファクタ（4 型構成 `SCEBasis` / `SCEDataset` /
`SCEFit` / `SCEModel`）の AtomsBase 統合を実装した際、adapter
(`src/AtomsBaseAdapter.jl`) は

```
AbstractSystem → TOML 形状の dict → Config4System
```

という経路を取った。kwarg コンストラクタ
`SCEBasis(; lattice, kd, ...)` も同様に `Config4System` を経由する。

`Config4System` は元々 TOML パース専用の中間型で、validation と
OffsetArray 構築のロジックを持つ。新 API ではすべての入力経路が
これを単一の parse / validation ターゲットとして共有する形になっている。

## 検討事項

理想の実装にとって `Config4System` は必要か、不要か。

## 現状の姿（2026-05-16 時点）

| 入力経路 | 経路 | 該当コード |
|---|---|---|
| TOML file | `TOML.parse → dict → Config4System` | `Magesty.jl` `SCEBasis(toml_path)` |
| Dict | `dict → Config4System` 直 | `Magesty.jl` `SCEBasis(input_dict)` |
| AtomsBase | `AbstractSystem → TOML-shape dict → Config4System` | `atomsbase_adapter.jl` `system_to_input_dict` (~30 行) + `_interaction_section` (~40 行) |
| kwargs | `kwargs → TOML-shape dict → Config4System` | `atomsbase_adapter.jl` `kwargs_to_input_dict` (~40 行) |

- `Config4System` は `ConfigParser` から export されているが、
  `Magesty` の公開 export list には含まれていない（半内部型）。
- Validation はパーサと癒着している：`parse_interaction_body1` 内で
  species 整合性、`parse_interaction_bodyn` 内で cutoff pair 完全性、
  `validate_system_parameters` でパラメータ境界をまとめて検査。
- Downstream (`Structure` / `Symmetry` / `Cluster` / `SALCBasis`) は
  `config::Config4System` を受け取り、フィールドを読むだけで再検証しない。
- `Config4System` は god config：system / symmetry / interaction /
  structure の 4 領域を 1 struct に詰めている。

## 案 A: Config4System を単一ゲートとして残す（pragmatic）

現状維持。`atomsbase_adapter.jl` の dict ビルダーも維持。
docstring と「internal-only」の明記でスコープを伝える。

**Pros**
- Validation は既に集約されているため再実装不要。
- 4 経路すべて動作・テスト済み。
- 1 コミット規模で完結。

**Cons**
- AtomsBase / kwargs adapter が「typed → dict → typed」の round-trip を
  踏んでおり、~150 行が dict シリアライザに費やされる。
- 5 番目の入力経路（CIF, JSON, etc.）を追加する際にも同じ dict
  ビルダーパターンを踏襲する負担が残る。
- god config が解消されない（`Symmetries.jl` などが必要のない
  `lattice_vectors` まで含む `Config4System` を持ち回る）。

## 案 B: typed value gate へ移行（ideal、本筋）

`Config4System` を「dict 入力専用パーサ」に縮退させ、各入力経路が
typed な値オブジェクトを直接組み立てる方向。god config も同時に分解。

```
                       ┌──→ SystemSpec      (name, atoms, lattice, positions, kinds)
[TOML file] → parser ──┼──→ InteractionSpec (nbody, lmax, lsum, cutoffs)
[Dict]      → parser ──┤
[AtomsBase] → direct  ─┼──→ SymmetryOptions (tolerance, isotropy)
[kwargs]    → direct  ─┘
                          ↓ (each: validate at construction)
                       SCEBasis(SystemSpec, InteractionSpec, SymmetryOptions)
```

**変更点**
- 3 つの typed value (`SystemSpec`, `InteractionSpec`, `SymmetryOptions`)
  を導入。各 inner constructor が自分の領域の validation を実施
  （species 整合性は `InteractionSpec` 内で `kd` を引数に取る形に）。
- TOML / Dict 専用の `parse_toml_inputs(dict) -> (SystemSpec,
  InteractionSpec, SymmetryOptions)` を分離。dict を扱う唯一の場所。
- AtomsBase adapter は dict を経由せず `SystemSpec(; lattice = ...,
  positions = ...)` を typed に直接呼ぶ。`system_to_input_dict` /
  `kwargs_to_input_dict` / `_interaction_section` は削除（~150 行 →
  ~50 行）。
- Downstream は必要な値のみを引数に取る：
  `Symmetry(structure, opts::SymmetryOptions)`,
  `Cluster(structure, inter::InteractionSpec)` 等。

**Pros**
- 「typed → dict → typed」round-trip が消え、AtomsBase / kwargs の
  実装が直接的になる。
- Validation が領域ごとに分離され、パーサと値オブジェクトの責任が
  明確になる。
- 新しい入力経路を追加する際は typed value を直接構築するだけ。
- god config が解消され、各サブモジュールが必要最小限のデータを受け取る。
- DFTK.jl / AtomsBase.jl エコシステムの設計慣習と整合。

**Cons**
- 数日〜1 週間規模の中規模リファクタ。spec 化 + `refactor/<slug>`
  ブランチが必要。
- `Config4System` を参照しているテスト・tools (`tools/personal/*` 含む)
  の書き直しが必要。
- 内部型ではあるが API 表面の破壊的変更（`Magesty.Config4System` を
  名前で参照しているコードは要修正）。

## 案 B の深掘り（2026-05-16 追記）

### B-1. typed value の定義案

| 型 | フィールド | validation 責任 |
|---|---|---|
| `SystemSpec` | `name::String`, `num_atoms::Int`, `kd_name::Vector{String}`, `kd_int_list::Vector{Int}`, `lattice_vectors::Matrix{Float64}` (3×3), `x_fractional::Matrix{Float64}` (3×N), `is_periodic::Bool` | lattice 行列形状、num_atoms と vector 長一致、kd_int_list の値域 ⊂ 1:length(kd_name)、x_fractional 範囲（periodic なら [0,1)） |
| `InteractionSpec` | `nbody::Int`, `body1_lmax::OffsetVector{Int}` (species ごと), `bodyn_lsum::OffsetArray{Int,2}` (body × species 等)、`bodyn_cutoff::OffsetArray{Float64,3}` (body × species × species) | nbody ≥ 1、body1_lmax のキーが kd_name と完全一致、bodyn_cutoff のペアが完全・対称、lsum 非負 |
| `SymmetryOptions` | `tolerance_sym::Float64`, `isotropy::Bool` | tolerance > 0 |

`InteractionSpec` の inner constructor は `kd_name::Vector{String}` を引数で受け取り
species 整合性検査を行うが、フィールドとしては保持しない（型レベルでは SystemSpec
に依存しない）。

### B-2. dict パーサの位置づけ

新ファイル `src/InputSpecs.jl`（仮）に 3 つの spec 型を集約。TOML/dict 経路専用の

```julia
parse_toml_inputs(dict::AbstractDict) -> (SystemSpec, InteractionSpec, SymmetryOptions)
```

を提供。これが dict を扱う**唯一の場所**で、ワイルドカード展開もここで完結する。
既存の `parse_interaction_body1` / `parse_interaction_bodyn` / `validate_system_parameters`
のロジックは各 spec の inner constructor に分散移植する（廃棄ではなく再配置）。

### B-3. 入力経路 4 種の再実装方針

| 経路 | 新実装 | 削減対象 |
|---|---|---|
| TOML | `TOML.parsefile → parse_toml_inputs(dict)` | なし（薄いまま） |
| Dict | `parse_toml_inputs(dict)` 直 | なし |
| AtomsBase | `SystemSpec(; ...)` + `InteractionSpec(; kd_name, ...)` を直接構築 | `system_to_input_dict` (~30 行) + `_interaction_section` の AtomsBase 分岐 (~20 行) |
| kwargs | 同様に typed を直接構築 | `kwargs_to_input_dict` (~40 行) + `_interaction_section` の kwargs 分岐 (~20 行) |

合計 ~110 行の dict ビルダー削減 + ~50 行の typed コンストラクタ追加 → 正味 ~60 行減。

AtomsBase / kwargs 経路でワイルドカードを使えるよう、共有 helper
`expand_wildcards(kd_name, entries) -> Dict` を提供（実体は `parse_toml_inputs` の
内部関数）。

### B-4. Downstream リファクタの順序

god config 分解は段階的に進められる：

1. `SystemSpec` / `InteractionSpec` / `SymmetryOptions` を追加（`Config4System` と並存させる）。
2. `parse_toml_inputs` を実装し、`Config4System` の構築経路を typed 経由に切り替え。
3. AtomsBase / kwargs adapter を typed 直構築に切り替え（dict ビルダー削除）。
4. Downstream を 1 つずつ移行: `Structure` → `Symmetry` → `Cluster` → `SALCBasis` → `Fitting` / `XMLIO`。
   各ステップで `Config4System` のフィールド参照を spec 参照に置換。
5. `Config4System` の参照が消えたら型ごと削除。

各ステップが独立コミットになり得る。

### B-5. テスト・tools 影響

- `test/component_test/test_ConfigParser.jl`: `Config4System` 直接構築の 8 箇所を
  `parse_toml_inputs(dict)` または各 spec の直接構築に置き換え（または `test_input_specs.jl`
  に統合）。
- `test/benchmark_salcbasis_hotspots.jl` / `test/develop_tmp/`: 同様に書き換え（develop_tmp は CI 外）。
- `tools/personal/*`: 必要に応じて移行。

### B-6. ワイルドカード species 指定（新機能）

ALAMODE `&cutoff` の語法を参考に、`bodyn_cutoff` と `body1_lmax` にワイルドカード
を導入する。`bodyn_lsum` は対象外。

**構文**:
- `bodyn_cutoff`: `"Fe-Ni"`, `"Fe-*"`, `"*-Fe"`, `"*-*"` の 4 形態。
- `body1_lmax`: `"Fe"`, `"*"` の 2 形態。

**Override 規則（specificity 勝ち）**:

| 優先度 | pair 例 | 単一 species 例 |
|---|---|---|
| 高 | `Fe-Ni`（完全一致） | `Fe`（完全一致） |
| 中 | `Fe-*` または `*-Fe`（片側ワイルドカード） | — |
| 低 | `*-*`（全ワイルドカード） | `*`（全ワイルドカード） |

- 同優先度キーが衝突する場合（例: `Fe-*` と `*-Ni` が両方 `Fe-Ni` をカバーし値が異なる）
  はエラー。ユーザに明示指定を促す。
- pair は無順序: `Fe-Ni` と `Ni-Fe` は同義。dict パーサで正規化（species 名の辞書順）
  してから優先度判定する。同義キーの両指定はエラー。
- 順序依存（TOML dict の hash 順）を避けるため、ALAMODE のテキスト順ベースではなく
  specificity-based を採用。記法は ALAMODE 互換だが override 解決は本パッケージ独自。

**展開タイミング**: `parse_toml_inputs(dict)` 内で全 pair / 全 species を具体形に展開
してから `InteractionSpec` を構築する。`InteractionSpec` の不変条件は従来通り
「全 species・全ペアが具体的に埋まっている」を保つ。

**TOML スキーマ**: 既存の `cutoff."Fe-Fe" = 5.0` 形式を踏襲し、`"*-*"` 等のキーを
許容（後方互換）。例:

```toml
[interaction.body2]
cutoff."*-*"  = 8.0
cutoff."Mg-O" = 10.0  # specificity 勝ちで Mg-O のみ 10.0

[interaction.body1]
lmax."*"  = 2
lmax."Fe" = 4         # Fe のみ 4、他全 2
```

**Validation の変更点**:
- 「全ペア完全性」: 展開後に `(kd[i], kd[j])` for `i ≤ j` の全組み合わせが埋まっている
  ことを `InteractionSpec` inner constructor で検査（従来通り）。
- 「specificity 衝突」: 展開時に検出して `parse_toml_inputs` でエラー。
- 「未定義 species 名」: `Fe-Xx` のような未知 species を含むキーはエラー（`*` のみ例外）。

### B-7. クラスタ定義の明文化

ワイルドカード機能の説明上、`cutoff` の意味を user-facing ドキュメントに明記する
必要がある。**既存実装** (`src/Clusters.jl` の `is_within_cutoff` および
`generate_clusters`) は：

> n-body interaction cluster (n ≥ 2) は、構成原子の **全 2-combinations** で各ペアの
> 距離が当該 body のペア cutoff (`cutoff[body, kd_i, kd_j]`) 以下を満たすときに採用
> される。

という定義を採っている。これは ALAMODE `&cutoff` の定義（cubic 以上も含めすべて
「全ペアが各 cutoff 内」で判定）と同じ。コードは変更せず、`docs/src/input_keys.md`
に明記する。

### B-8. 残された未決事項（spec 化時に詰める）

- `SystemSpec.lattice_vectors` を `SMatrix{3,3,Float64}` にするか `Matrix{Float64}` の
  ままにするか（StaticArrays 化はホットパスでの効果次第）。
- `parse_toml_inputs` のエラー集約方法（現状の 1 箇所 throw vs 複数エラー集約）。
- `body1_lmax` / `bodyn_cutoff` 以外への将来のワイルドカード拡張（`bodyn_lsum` 等）
  は spec の Future work に記載しておく。

## 結論（2026-05-16）

**案 B に進む**。当初は「5 番目の入力経路追加時に再判断」と保留していたが、
ワイルドカード species 指定の新機能要請により、dict パース層に substantive な責務
（specificity-based override 解決）が加わることが判明。案 A での後付けは god config の
弊害を温存したまま dict ビルダーが肥大するため、案 B 移行と一度でやるほうが筋がよい。

spec フォルダ `docs/specs/260516-typed-input-spec/` で詳細設計を確定させてから
`refactor/typed-input-spec` ブランチで実装する。
