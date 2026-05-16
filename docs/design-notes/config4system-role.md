# Config4System の役割の見直し

**Status**: 結論保留（2026-05-14 提案 / 2026-05-16 整理）。
今は何もしない。入力経路が 5 つ目を必要とするタイミングで再判断する。

**関連**: [`docs/specs/260514-sce-public-api/`](../specs/260514-sce-public-api/)
（SCE 公開 API リファクタ。merged）

## 背景

SCE 公開 API リファクタ（4 型構成 `SCEBasis` / `SCEDataset` /
`SCEFit` / `SCEModel`）の AtomsBase 統合を実装した際、adapter
(`src/utils/atomsbase_adapter.jl`) は

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

## 結論（2026-05-16）

**理想は案 B**。現状の dict round-trip は「TOML-first だった頃の
副産物」であり、典型的な Julia エコシステム慣習からは外れている。
AtomsBase adapter の ~150 行 dict ビルダーが消え、god config も
解消される利点は大きい。

ただし **今すぐ着手はしない**：

- 4 入力経路すべて動作・テスト済みで、実利的な不便が出ていない。
- god config 分解は効果が大きいが波及範囲も大きく、他にも
  優先すべきタスクがある（performance backlog 等）。
- 5 番目の入力経路（CIF / JSON / database / ...）を追加する話が
  出たタイミングで「ついでに案 B に移行」が自然。

着手判断のトリガー：
- 新しい入力経路を 1 つでも追加する話が出たとき。
- god config の弊害（`Symmetry` / `Cluster` が無関係なフィールドを
  保持している点）が他の機能追加の妨げになったとき。
- AtomsBase adapter に新たな機能（kind サブラベル、Wyckoff 区別 等）
  を加える必要が出たとき — typed gate に移行した方が実装が楽になる。
