# リファクタリング候補スイープ R1-R11

**Status**: ほぼ完了（2026-05-13 開始、2026-05-14 完了）。R1-R7, R9, R10, R11 完了。R8 は Plan B (重複集約) 完了、Plan C (バグ修正) は `replace-sorted-container.md` の DataStructures.jl 置き換え案で吸収予定のため保留。

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

**Status**: **完了** (2026-05-14). branch なし → main に直接 fast-forward 予定。

**対象**: `src/utils/AngularMomentumCoupling.jl`

**実態確認結果**: 当初の見立て（「`build_all_complex_bases` / `build_all_real_bases` が両方を呼んでいる」）は誤り。実際の callers を grep した結果:

- `coeff_tensor_complex`: `build_all_complex_bases` (L389)、`build_all_real_bases` (L425)、テスト 6 箇所で使用。
- `coeff_tensor_complex_mindexed`: in-repo の caller ゼロ（src/ にも test/ にも tools/ にもなし）。**初期コミット `144e620` から存在するが当初から完全な dead code**。

**実装**: dead code 削除を採用（共通 DP コア抽出ではなく）。pre-release で内部利用ゼロ、テストゼロのため、削除は極めて低リスク。削除対象:
- `coeff_tensor_complex_mindexed` 関数本体（旧 L267–359）
- `using OffsetArrays` (旧 L28、本ファイル内で他に利用なし)
- `export` 一覧の `coeff_tensor_complex_mindexed,`

`OffsetArrays` パッケージ依存自体は `Clusters.jl` / `BasisSets.jl` / `ConfigParser.jl` が使うため Project.toml からは外せない。

**連動箇所**: なし（モジュール内に閉じている）。テストは既存の 428 件全パス。`make test-integration` も 155/155 全パス。

## 🟡 中優先度（保守性・スケーラビリティ）

### R6. `MySphericalHarmonics.jl` の buffered / non-buffered API の整理

**Status**: **完了** (Plan B, 2026-05-14). branch なし → main に直接 fast-forward 予定。

**対象**: `src/utils/MySphericalHarmonics.jl`

**実態確認結果**: buffered 関数は計 5 つ (`_legendre_pair_unsafe!`, `P̄ₗₘ`, `dP̄ₗₘ_unsafe`, `Zₗₘ_unsafe`, `∂ᵢZlm_unsafe`)。buffer サイズ要件のコメント / docstring が 5 箇所に散在しており、しかも `dP̄ₗₘ_unsafe` だけ要件が異なる (`l - |m|` vs 他 4 つの `l - |m| + 1`)。実態として `∂ᵢZlm_unsafe` が両者を 1 つのバッファで共有するため、caller は `l - |m| + 1` を確保する必要があり、`dP̄ₗₘ_unsafe` 単独要件の `l - |m|` は混乱の元。

**実装** (数値結果不変、docstring + boundscheck のみ):
1. モジュール docstring に "Buffer requirements" セクションを追加し、全 buffered 関数の要件を一箇所に集約。`l - |m| + 1` を統一基準として明記。`max_l + 1` をスレッドごとに確保する推奨パターンを記載。
2. private helper `_required_buf_size(l, m) = l - abs(m) + 1` を追加。
3. 5 つの buffered 関数に `@boundscheck checkbounds(buf, _required_buf_size(l, m))` 追加。`@inbounds` 環境では elided されるため hot path 性能は不変。
4. 散在していた buffer サイズコメント (`P̄ₗₘ` L132, `dP̄ₗₘ_unsafe` L154, `Zₗₘ_unsafe` docstring L482-490, `∂ᵢZlm_unsafe` docstring L738-745) を「See the module docstring's Buffer requirements section.」への参照に簡潔化。

**連動箇所**: 球面調和関数の規格化・符号は変更なし → `SphericalHarmonicsTransforms.jl` / `BasisSets.jl` / `test_sphericart_agreement.jl` への影響なし。テスト: `test-unit` 6203/6203、`test-integration` 155/155、`test-jet`、`test-aqua` 全パス。

### R7. `xml_io.jl` の XML タグ・属性のマジックストリング

**Status**: **完了** (commit `cc53c8c`, 2026-05-13). branch なし → main に直接 fast-forward。

**対象**: `src/utils/xml_io.jl`

`"Magesty"`, `"System"`, `"SALC"`, `"basis"` などのノード名・属性名、`"%.15e"` のフォーマット文字列が read/write 両側で hardcoded（連動箇所: SCE 係数の入出力で整合性が求められている部分）。

**実装**: モジュール冒頭に `TAG_*` / `ATTR_*` 定数ブロックを追加し、write/read 両方を経由させる形に統一。`@sprintf` はリテラルフォーマットしか受け付けないため、`fmt_lattice` / `fmt_fractional` / `fmt_tensor` の 1 行 wrapper 関数で間接化。`<$(TAG_*)>` interpolation を例外メッセージにも適用してスキーマ名変更時の整合性を自動化。XML 出力は byte-identical（integration tests で write/read round-trip 確認）。

**連動箇所**: integration tests (dimer, fept, fege など) が write_xml + build_sce_basis_from_xml を経由しているため自動カバー。専用 round-trip テストは未追加（不要と判断）。

### R8. `SortedContainer.jl` の 3 コンテナの実装重複

**Status**: **完了** (Plan B, 2026-05-14). branch なし → main に直接 fast-forward 予定。Plan C (バグ修正・push!/delete! 共通化) は別 PR 扱い。

**対象**: `src/common/SortedContainer.jl`

**実態確認結果**: `AbstractSortedVector` が抽象タグとして定義されていたが、generic メソッドは一切定義されておらず、`getindex` / `length` / `size` / `iterate` / `isempty` / `append!` / `findfirst` / `in` を 3 コンテナそれぞれが個別実装していた。`SortedVector` と `SortedUniqueVector` の `findfirst` / `in` は完全に同一の `searchsortedfirst` ベース実装。`deleteat!` / `clear!` も 2 コンテナで完全同一だが、`SortedCountingUniqueVector` には未定義 (`counts` Dict との同期が必要なため)。

**実装** (Plan B: pure refactor、振る舞い変更なし):
1. `AbstractSortedVector` の docstring を更新し、`data` フィールドの contract を明文化 (sorted な indexable コレクション、`Vector{T}` または別の `AbstractSortedVector{T}`)。
2. 以下を generic fallback として `AbstractSortedVector` に集約:
   - `getindex(asv, i::Int)` / `length(asv)` / `size(asv)` / `iterate(asv)` / `iterate(asv, state)` / `isempty(asv)`
   - `append!(asv, vec)` (`push!` を呼ぶ標準パターン)
   - `findfirst(asv, value)` / `in(value, asv)` (両者とも `searchsortedfirst` ベース)
3. 3 コンテナから上記の重複メソッドを削除。
4. `SortedCountingUniqueVector::in` のみ `counts` Dict ベース (O(1)) の固有実装を override として保持。
5. `deleteat!` / `clear!` は **意図的に集約しない**: `SortedCountingUniqueVector` で適用すると `counts` Dict と `data` が不整合になるため、`SortedVector` / `SortedUniqueVector` 個別実装のまま残す。これらの正しい統合は Plan C で扱う。

**ハマりポイント**: `in(val, scv)` の override シグネチャで `val` の型を書かなかったため generic `in(value::T, asv::AbstractSortedVector{T})` と method ambiguity が発生。`val::T` を明示することで解消。

**残課題 (Plan C で別 PR)**:
- `SortedCountingUniqueVector::clear!` / `deleteat!` 未定義 (現状は MethodError、`data` も `counts` も両方クリアする実装が必要)。
- `SortedVector::delete!` / `SortedUniqueVector::delete!` が sorted 性を活かさず O(N) (`searchsortedfirst` で O(log N) 可能)。
- `SortedUniqueVector::copy` が無駄な再 sort + unique を実行。
- `push!` の sorted insertion 検索部分の共通化。

**連動箇所**: なし（モジュール内に閉じている、外部 API シグネチャ不変）。テスト: `test-unit` 6203/6203、`test-integration` 155/155、`test-jet`、`test-aqua` 全パス。差分: 386 → 356 行 (-30 行)。

### R9. 公開 API の docstring 不揃い

**Status**: **完了** (2026-05-14). branch なし → main に直接 fast-forward 予定。

**対象**: `src/Magesty.jl` の export 一覧および周辺 getter。

**実態確認結果**: `src/Magesty.jl` の export を起点に audit。多くの docstring は既に Julia 標準形式（# Arguments / # Returns / # Examples）に揃っていたが、以下が不揃いだった:

- `get_j0` / `get_jphi` / `get_j0_jphi` (qualified access のみ、`Magesty.get_j0(sc)`): 1 行サマリのみだった
- `build_sce_basis`: 1 行サマリのみだった
- `SpinCluster(system, input_dict, spinconfig_list)` 3 引数 form: 2 行のみだった
- `predict_energy`: # Examples 欠落
- `install_tools`: 散文形式、# Arguments / # Returns 欠落
- `fit_sce_model` の docstring 例: `weight=0.7` と書かれていたが `weight` は positional 引数のため Julia ではエラーになる呼び出し（バグ）

**実装**: 上記すべてを Julia 標準形式に統一。`fit_sce_model` の例は `0.7`（位置引数）に修正。`get_j0` 系の signature 行は `get_j0(sc::SpinCluster) -> Float64` のように戻り値型まで明示。

**連動箇所**: なし（docstring のみ）。テスト: `test-unit` 6203/6203、`test-jet`、`test-aqua` 全パス。

### R10. `Basis.jl` `CoupledBasis` constructor の不変条件検証の散在

**Status**: **完了** (2026-05-14). branch なし → main に直接 fast-forward 予定。

**対象**: `src/types/Basis.jl`

**実態確認結果**:
- 当初想定の「docstring で atoms が sorted と書かれている」は **事実誤認**。struct docstring に sorted 文言なし、inner constructor (旧 L75) に `# Check if atoms are sorted` という dead comment があるのみ。caller の `tesseral_coupled_bases_from_tesseral_bases` も atoms を sort せずに渡しており、不変条件として強制すると既存挙動が壊れる → atoms sort 検証は **追加しない**。
- struct docstring 2 箇所 (`CoupledBasis`, `CoupledBasis_with_coefficient`) で `Lseq` の長さを "length N-1" と記述していたが、inner constructor の実装は `N-2` をチェック → docstring の誤記。N-2 が正 (`Lf` は別フィールド)。
- `CoupledBasis_with_coefficient` の inner constructor は `length(Lseq)` と `length(atoms)` のチェックが完全に欠落していた（`CoupledBasis` 側にはある）。
- `CoupledBasis` の N=1 early return が `length(atoms) == 1` と `ndims(coeff_tensor) == 2` のチェックをスキップしていた。

**実装**:
1. docstring の "length N-1" → "length `max(0, N-2)`" に修正 (両 struct)。
2. `CoupledBasis_with_coefficient` constructor に `length(Lseq) == max(0, N-2)` と `length(atoms) == N` チェック追加。
3. `CoupledBasis` の N=1 early return を廃止し、共通パスでバリデーション通過。N=1 のときも `max(0, N-2) = 0`、`length(atoms) == 1`、`ndims(coeff_tensor) == 2` を強制。
4. エラーメッセージフォーマットを `"<field> must be <constraint> = <expected>; got <actual> (N=<N>)"` に統一。
5. dead comment `# Check if atoms are sorted` を削除。

**連動箇所**: なし（モジュール内に閉じている、外部 API シグネチャ不変）。テスト: `test-unit` 6203/6203、`test-integration` 155/155、`test-jet`、`test-aqua` 全パス。

### R11. `EnergyTorque.jl` の `(4π)^(n_C/2)` スケーリングの重複と物理意味未注釈

**Status**: **完了** (2026-05-14). branch なし → main に直接 fast-forward 予定。

**対象**: `src/utils/EnergyTorque.jl` `calc_energy` / `calc_torque`, `src/Optimize.jl` `build_design_matrix_energy` / `build_design_matrix_torque` / `predict_energy`

**実態確認結果**: `(4π)^(n_C/2)` の式が 5 箇所に散在 (`Optimize.jl` の L363 / L514 / L887、`EnergyTorque.jl` の L59 / L134)。リテラルも `4*pi` と `4π` が混在。物理由来 (技術ノート参照) の docstring はどこにも書かれていない。SCE 基底関数の `1/√(4π)` 正規化を打ち消すための係数であることが caller では一切わからない状態。

**参考**: `SpinClusterMC.jl` (`src/JPhiMagestyCarlo.jl:141-144`) では既に `@inline _cluster_scaling(n_sites::Integer)::Float64 = (4 * pi)^(n_sites / 2)` として private helper に集約済み。同じパターンを Magesty 側にも適用。

**実装** (Plan B: helper 抽出のみ、式不変):
1. `Optimize.jl` の冒頭付近に private helper を追加:
   ```julia
   @inline _cluster_scaling(n_sites::Integer)::Float64 = (4π)^(n_sites / 2)
   ```
   docstring に物理由来 (技術ノート参照) を明記。export しない。
2. 5 箇所すべてを `_cluster_scaling(n_C)` 呼び出しに置換。`EnergyTorque.jl` 側は `Optimize._cluster_scaling(n_C)` で参照 (既存の `using ..Optimize` で十分)。
3. 散在していた `# (√(4π))^{n_C}` コメントは helper の docstring に集約されるため削除。

**数値結果**: 完全に identical (式不変、Float64 で `4*pi` も `4π` も同じ値)。integration tests (energy / torque を baseline と比較) で確認。

**スコープ外 (将来の整理候補)**:
- `EnergyTorque.calc_energy` ≈ `Optimize.predict_energy`、`EnergyTorque.calc_torque` と `Optimize.build_design_matrix_torque` の構造重複は本 R11 では触らない。`sce-public-api.md` の predict 系 API 整理と整合させて対応 (`EnergyTorque.jl` を deprecate or 統合)。

**連動箇所**: 物理規約は不変。CLAUDE.md「連動箇所」セクションに「`_cluster_scaling` の規約は技術ノートと整合」項目を追加する候補だが、helper 化により drift リスクは大幅減 (5 箇所 → 1 箇所)。テスト: `test-unit` 6203/6203、`test-integration` 155/155、`test-jet`、`test-aqua` 全パス。

## 🟢 低優先度（参照のみ）

- **`BasisSets.jl` legacy SHProduct path 群**: ~~`construct_basislist`（L957–1000 付近）~~ **完了** (Phase 1 + Phase 2, 2026-05-14). Phase 1: `construct_basislist` を起点とする 9 関数 (`push_unique_body!`, `listup_basislist`, `map_atom_l_list`, `classify_basislist`, `is_translationally_equiv_basis`, `projection_matrix` (legacy), `proj_matrix_a_symop`, `operate_symop`, `corresponding_idx`) を `BasisSets.jl` から削除 (~500 行)。Phase 2: 連動して `AtomicIndices` モジュール全体が dead になったので、`src/types/AtomicIndices.jl` (217 行) と `test/component_test/test_AtomicIndices.jl` (249 行) を削除。`SHSiteIndex` / `SHProduct` / `LinearCombo` / `inner_product` / `replace_atom` / `product_shsiteindex` / `get_atom_l_list` / `shsiteindex_singleatom` がすべて消滅。BasisSets / Optimize / Magesty / runtests の `using ..AtomicIndices` / `include` も同時に整理。
- **`SpinConfigs.show`（L175–197）と `Clusters.print_cluster_stdout`（L506–566）**: 出力フォーマットの統一余地。`PrettyTables.jl` 検討の余地もあるが新規依存追加なので保留。
- **`AtomCells.jl` `AtomCell` と `AtomicIndices.jl` `SHSiteIndex`**: `isless` / `==` / `hash` の実装パターンが重複。**注意**: 物理意味が違うので無理に共通化するとバグの温床。慎重に。
- **`common/version.jl`**: ~~hardcoded version 文字列。Project.toml から動的取得する仕組みの検討余地。~~ **完了** (commit `95269de`, 2026-05-13). `pkgversion(@__MODULE__)` で Project.toml から動的取得する形に変更。bump 時の編集箇所は Project.toml の 1 行のみ。`test/component_test/test_Version.jl` で drift を検証。
- **`Magesty.jl` `print_header`（L743–757）**: version / Julia info / threads / timestamp の 4 責務。verbose ロギングを整備する際に分割。

## Optimize.jl 関連項目（クロスリファレンス）

`alpha` パラメータの "currently unused"、`_fit_sce_model_internal` の `isa` 分岐、`elastic_net_regression` の責務混在、`Optimizer` コンストラクタの曖昧な API は、`estimator-dispatch.md`（完了済み）で詳細カバー済み。重複記述しない。
