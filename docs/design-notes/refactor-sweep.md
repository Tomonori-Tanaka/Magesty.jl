# リファクタリング候補スイープ R1-R11

**Status**: 進行中（2026-05-13 開始、R1 / R2 / R3 / R4 / R7 完了）

**目的**: `src/` 全体（~9.4k 行）を対象に「構造的に整理すべき」候補を洗い出した結果。Optimize.jl は `estimator-dispatch.md` で個別カバー済みなので、本ファイルでは触れない。

**着手判断**: 個別の項目を実装する前に CLAUDE.md の方針に従い、中規模リファクタリングは spec 化（`docs/specs/[YYMMDD]-[slug]/`）が必要。バグ修正やドキュメント追加で完結する小規模項目は spec なしで進めてよい。「要 verify」マーク付きの項目は、着手前に該当ファイルを読んで実態を確認する。

## 🔴 高優先度（API・数値結果・既知バグに関わるもの）

### R1. `System` / `build_sce_basis` / `SpinCluster` の構築シーケンス重複

**Status**: **完了** (commit `2ffc2c9`, 2026-05-13). branch なし → main に直接 fast-forward。

**対象**: `src/Magesty.jl`

`System(input_dict)`（L132–148）、`build_sce_basis`（L172–179）、`SpinCluster(input_dict)`（L286–307）の 3 箇所で、`Config4System → Structure → Symmetry → Cluster → BasisSet` の同じ 5 ステップが繰り返されていた。さらに `build_sce_basis_from_xml`（L205–230）は BasisSet 部分だけ XML 読み込みに置換した 4 つ目のバリアント。

**実装**: private helper `_build_structure_skeleton(config; verbosity)` を抽出し `(structure, symmetry, cluster)` を返す形に。BasisSet は caller が 1 行で追記する（XML パスは置換しやすい）。当初案の `basisset_override` kwarg は採用せず、caller 側で分岐する方が `build_sce_basis_from_xml` の verbose 出力順序が保たれる。

**連動箇所**: なし（モジュール内に閉じている）。

### R2. `write_xml` の 3 シグネチャと位置/キーワード引数の混在

**Status**: **完了** (commit `d9ae9ba`, 2026-05-13). Spec: `docs/specs/260513-write-xml-api/`. Breaking change: 4-arg public method 削除。

**対象**: `src/Magesty.jl` `write_xml`（L477–540 付近）

調査で判明した実状: 公開シグネチャは 3 つではなく **4 つ** 存在していた（`SpinCluster`, `System`, `structure+symmetry+basis_set+optimize` の 3 つに加え、内部 `XMLIO.write_xml` の 4-arg helper）。`write_xml(system, "out.xml", false)` のような呼び出しで第 3 引数の意味が型ディスパッチによって変わり、ユーザーから読み取れない問題があった。

**実装**: 当初案の「完全 kwarg 化 (`write_xml(target; filename, write_jphi=true)`)」は Julia エコシステム慣習 (`Base.write`, `CSV.write`, `JLD2.save` 全て filename 位置引数) と整合しないため採用せず。代わりに公開 4-arg form を削除し、2 メソッドに統一:

```julia
write_xml(system::System, filename = "jphi.xml")
write_xml(sc::SpinCluster, filename = "jphi.xml"; write_jphi = true)
```

XMLIO 側の内部 helper は据置（`SpinCluster` overload の実装 backend）。In-repo の唯一の 4-arg caller (`docs/src/examples.md`) は `SpinCluster(structure, symmetry, cluster, basisset, optimizer)` でラップしてから `write_xml(sclus, ...)` を呼ぶ形に rewrite。

数値: XML byte-identical (`dimer.xml`, `fept/scecoeffs.xml` で baseline 比較済み)。

**連動箇所**: `xml_io.jl` の内部 helper、`docs/src/examples.md` (L186)。

### R3. `BasisSets.jl:554` の TODO とコメントアウト処理（既知バグ）

**Status**: **完了** (commit `833b153`, 2026-05-13). branch なし → main に直接 fast-forward。

**対象**: `src/BasisSets.jl` `listup_coupled_basislist` 周辺

L554 に「Remove this exceptional handling when the bug is fixed in the projection matrix construction」のコメントとともに、コメントアウトされた `if sort(l_vec) == [1, 3]; continue; end` が残存していた。Bug guard は `8b389b4` (Feb 2026) で追加、`3fc4180` (Mar 2026) の `find_translation_atoms` refactor でコメントアウトに変更（"presumed fixed"）。回帰テストは未追加だった。

**実装**: 2-atom BCC + `body2.lsum = 4` で `l_vec ∈ {[1,3], [2,2], [3,1]}` を経由する回帰テストを `test/component_test/test_BasisSets_l13_regression.jl` に新設し、`projection_matrix_coupled_basis` まで含めた full BasisSet build が通ることを確認。その後コメントアウトコードと TODO を削除。

**連動箇所**: テスト追加済み（`test_BasisSets_l13_regression.jl`、4 テスト）。

### R4. `ConfigParser.jl` のデフォルト値・検証ロジックの 3 箇所散在

**Status**: **完了** (commit `73edde3`, 2026-05-13). branch なし → main に直接 fast-forward。

**対象**: `src/utils/ConfigParser.jl`

**実装上の発見**: `VALIDATION_RULES_SYSTEM` は定義されていたが `Config4System` の constructor が validation 関数を一切呼んでおらず、**完全な dead code** だった。重複していたのは「`if nbody < 1; throw(...)` を imperative に書く」のと「`ValidationRule(:nbody, x -> x > 0)`」の 2 表現があった点。

**実装**: 当初案（"declarative テーブルから注入"）ではなく、Julia の慣習に合わせた逆方向の整理を採用:
1. `ValidationRule` 構造体・`VALIDATION_RULES_*` 定数・`_apply_validation_rules` を **削除**（型消去された validator フィールド・cross-field 検証不可・適用先 2 箇所しかない、と理由は重なる）。
2. `validate_system_parameters` / `validate_optimize_parameters` を **imperative な inline check** で書き直し（DFTK.jl, AtomsBase.jl 等の標準パターン）。
3. 必須セクションチェックだけは `_check_required_sections` + `REQUIRED_SECTIONS_*` 定数で共通化（Config4System / Config4Optimize で再利用）。
4. 副次的挙動変更: `name = ""` のような無効入力が以前は dead validation で通っていたのが、今は意図通り `ArgumentError` を throw。テスト 2 件追加 (`name=""`, `nbody=0`)。

**連動箇所**: `Config4System`, `Config4Optimize`、TOML スキーマドキュメント（変更なし）。

### R5. `AngularMomentumCoupling.jl` の DP ロジック重複

**対象**: `src/utils/AngularMomentumCoupling.jl` `coeff_tensor_complex`（L137–212）と `coeff_tensor_complex_mindexed`（L279–359 付近、要 verify）

3 ステージの CG 漸化 DP（N ≥ 3）がほぼ同じアルゴリズムで 2 度実装されている。indexing 方式（`CartesianIndices` vs m-direct `OffsetArray`）のみ違う。

**改善案**: 共通の DP コアを 1 関数に抽出し、indexing 方式を wrapper で選択する。

**連動箇所**: `build_all_complex_bases`, `build_all_real_bases` が両方を呼んでいる。

## 🟡 中優先度（保守性・スケーラビリティ）

### R6. `MySphericalHarmonics.jl` の buffered / non-buffered API の整理

**対象**: `src/utils/MySphericalHarmonics.jl`

`Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` で buffered 版と non-buffered 版が並列に存在。buffer サイズ要件（`length(buf) >= l - |m| + 1` 等）の docstring が散在しており、新規呼び出し時に必要サイズの根拠が辿りにくい。

**改善案**: buffer サイズ仕様を 1 箇所（例: モジュール冒頭の docstring または `BUFFER_SIZE_*` 定数）に集約し、各 buffered 関数の docstring から参照する。可能なら `@boundscheck` で buffer サイズの assertion を追加。

### R7. `xml_io.jl` の XML タグ・属性のマジックストリング

**Status**: **完了** (commit `cc53c8c`, 2026-05-13). branch なし → main に直接 fast-forward。

**対象**: `src/utils/xml_io.jl`

`"Magesty"`, `"System"`, `"SALC"`, `"basis"` などのノード名・属性名、`"%.15e"` のフォーマット文字列が read/write 両側で hardcoded（連動箇所: SCE 係数の入出力で整合性が求められている部分）。

**実装**: モジュール冒頭に `TAG_*` / `ATTR_*` 定数ブロックを追加し、write/read 両方を経由させる形に統一。`@sprintf` はリテラルフォーマットしか受け付けないため、`fmt_lattice` / `fmt_fractional` / `fmt_tensor` の 1 行 wrapper 関数で間接化。`<$(TAG_*)>` interpolation を例外メッセージにも適用してスキーマ名変更時の整合性を自動化。XML 出力は byte-identical（integration tests で write/read round-trip 確認）。

**連動箇所**: integration tests (dimer, fept, fege など) が write_xml + build_sce_basis_from_xml を経由しているため自動カバー。専用 round-trip テストは未追加（不要と判断）。

### R8. `SortedContainer.jl` の 3 コンテナの実装重複

**対象**: `src/common/SortedContainer.jl`

`SortedVector`, `SortedUniqueVector`, `SortedCountingUniqueVector` が `AbstractSortedVector` を継承する設計だが、`push!` / `append!` / `delete!` の本体ロジックが重複している（要 verify）。

**改善案**: 共通の sorted insertion / deletion 戦略を base に統合し、uniqueness / counting を直交した拡張として実装する。

### R9. 公開 API の docstring 不揃い

**対象**: 主に `src/SpinConfigs.jl` `get_j0`, `get_jphi`, `get_j0_jphi`（L683–701 付近）と export 一覧全体

CLAUDE.md は「エクスポートされる API には明示的な型アノテーションと docstring を付ける」「Julia 標準形式（`# Arguments` / `# Returns` / `# Examples`）」を要求しているが、一部の getter で Arguments / Returns セクションが欠落。

**改善案**: `src/Magesty.jl` の export 一覧（L77–84）を起点に audit し、全 export 関数の docstring を Julia 標準形式に統一する。

### R10. `Basis.jl` `CoupledBasis` constructor の不変条件検証の散在

**対象**: `src/types/Basis.jl:43–79`

docstring では `atoms` が "sorted" であるべきと書かれているが、inner constructor で assertion されていない（要 verify）。N=1 special case の handling、Lseq 長チェック、エラーメッセージのフォーマットが不揃い。

**改善案**: すべての不変条件を inner constructor 1 箇所に集約し、エラーメッセージのフォーマットを統一する。

### R11. `EnergyTorque.jl` の `(4π)^(n_C/2)` スケーリングの重複と物理意味未注釈

**対象**: `src/utils/EnergyTorque.jl` `calc_energy`（L59 付近）, `calc_torque`（L134 付近）, および `src/Optimize.jl` `build_design_matrix_energy` / `build_design_matrix_torque` 内の同様の式

物理由来（technical_notes 参照）の docstring が欠落しているうえ、3 箇所以上に同じ式が散らばっている。スケーリングの規約が変わったら全箇所同期が必要だが、現状は静かに drift しうる。

**改善案**: 定数または `_sce_scaling_factor(n_C)` helper として 1 箇所に集約し、由来（[Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/)）を docstring に明記。

**連動箇所**: CLAUDE.md「連動箇所」セクションに追加する候補。

## 🟢 低優先度（参照のみ）

- **`BasisSets.jl` `construct_basislist`（L957–1000 付近）**: "backward compatibility" コメント付きで、内部呼び出しがあるかどうか要 verify。なければ削除。
- **`SpinConfigs.show`（L175–197）と `Clusters.print_cluster_stdout`（L506–566）**: 出力フォーマットの統一余地。`PrettyTables.jl` 検討の余地もあるが新規依存追加なので保留。
- **`AtomCells.jl` `AtomCell` と `AtomicIndices.jl` `SHSiteIndex`**: `isless` / `==` / `hash` の実装パターンが重複。**注意**: 物理意味が違うので無理に共通化するとバグの温床。慎重に。
- **`common/version.jl`**: ~~hardcoded version 文字列。Project.toml から動的取得する仕組みの検討余地。~~ **完了** (commit `95269de`, 2026-05-13). `pkgversion(@__MODULE__)` で Project.toml から動的取得する形に変更。bump 時の編集箇所は Project.toml の 1 行のみ。`test/component_test/test_Version.jl` で drift を検証。
- **`Magesty.jl` `print_header`（L743–757）**: version / Julia info / threads / timestamp の 4 責務。verbose ロギングを整備する際に分割。

## Optimize.jl 関連項目（クロスリファレンス）

`alpha` パラメータの "currently unused"、`_fit_sce_model_internal` の `isa` 分岐、`elastic_net_regression` の責務混在、`Optimizer` コンストラクタの曖昧な API は、`estimator-dispatch.md`（完了済み）で詳細カバー済み。重複記述しない。
