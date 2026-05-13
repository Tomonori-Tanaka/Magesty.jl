# Design Note

改善案・設計メモを記録するファイル。実装済みの履歴は `.claude/bench_log.md` と `git log` を参照。

---

## 調査結果: SpheriCart.jl 採用の可否（2026-05-11）

**結論: Magesty.jl では SpheriCart を採用しない**（SpinClusterMC とは逆の結論）。

### 背景
SpinClusterMC 側調査（同日付）で Magesty `Zₗₘ_unsafe` ≡ SpheriCart `Y_l^m_real`
が bit-exact 一致することは確認済み。Magesty 本体でも採用が有利か調査。

### 動作確認
- `make test-sphericart` → `SpheriCart vs MySphericalHarmonics agreement` で
  700 / 700 PASS（`test/component_test/test_sphericart_agreement.jl`、lmax=4 全 (l,m)、7 単位ベクトル）
- 数値同一性は完全（SpinClusterMC 側調査と同じ）

### 性能実測（`make bench-optimize-sphericart` 系）

`test/benchmark_optimize_sphericart.jl` は v1（drop-in: per-atom `compute`）と
v2（per-spinconfig 事前計算キャッシュ）両方を実装済み。median 比較：

**fege_2x2x2（64 atoms, lmax = 1, 146 SALCs, 100 spinconfigs）**

| 関数 | baseline (Magesty) | v1 (SpheriCart drop-in) | v2 (SpheriCart 事前計算) |
|---|---|---|---|
| `design_matrix_energy_element` | **54.0 μs** | 133 μs (2.46× slow) | 113 μs (2.09× slow) |
| `calc_∇ₑu` | **8.4 μs** | 23.5 μs (2.80× slow) | 21.5 μs (2.55× slow) |
| `build_design_matrix_torque`（1 config） | **208 ms** | — | 503 ms (2.42× slow) |

**fept_tetragonal_2x2x2（16 atoms, lmax = 2, 31 SALCs, 30 spinconfigs）**

| 関数 | baseline | v1 | v2 |
|---|---|---|---|
| `design_matrix_energy_element` | **27.7 μs** | 2.40× slow | 2.07× slow |
| `calc_∇ₑu` | **7.7 μs** | 3.11× slow | 2.71× slow |
| `build_design_matrix_torque` | **6.56 ms** | — | 2.39× slow |

メモリも v2 は torque 行列で **785 MiB** vs baseline **293 MiB**（2.7×）。

### なぜ SpinClusterMC とは逆の結果か

アクセスパターンの違い：

- **SpinClusterMC**: sweep 中の atom cache は「ある atom について全 (l, m) up to L_max」を
  欲しい。SpheriCart の SVector 一括返しが access pattern と完全一致 → **10-20× 高速化**。
- **Magesty `design_matrix_energy_element` / `calc_∇ₑu`**: ループは
  `for (site_idx, atom) in translated_atoms; l = cbc.ls[site_idx]; for m_idx in 1:2l+1`
  と、**site ごとに固定の `l`** で 2l+1 個の m 値しか要らない。SpheriCart は
  L_max まで全 (l, m) 計算するため、その他の (l', m') 計算は完全に無駄。
  L_max=2 でも (L+1)² = 9 個のうち欲しいのは ~1-3 個程度 → SIMD 並列の節約より
  「不要 (l, m) を計算する分のオーバーヘッド」が勝る。

v2（事前計算キャッシュ）でも勝てない理由：
- per-spinconfig の `Y_cache, dY_cache` 構築は cheap だが、各 `itrans × site` 反復で
  `sh_values = [[Y_cache[...] for mi in 1:2l+1] for si in 1:N]` と **inner Vector を
  毎回ヒープ確保**する構造になっており、確保コストが Magesty 内蔵 `Zₗₘ_unsafe` の
  値計算コストを上回る。

### Magesty 側で SpheriCart を活かす条件（参考）

理論上 SpheriCart 有利になる構造は：
1. すべての SALC の `cbc.ls` 最大値 ≈ 個別 `l` 値の平均（i.e., 全 (l,m) を毎回使う）
2. inner loop の `Vector{Float64}` 確保を完全に排除（`MMatrix` / 線形 indexing）
3. spinconfig 単位の cache を `Matrix{Float64}` のまま row として読む

しかし現状 Magesty の SALC は cluster ごとに `ls` が混在（1〜lmax まで様々）、
高 lmax SALC は SALC 数全体の小さな割合。restructure コストに見合う改善は望めない。

### 採用条件下での将来再評価指針

以下が同時に成立する場合は再検証する価値あり：
- L_max ≥ 4 の SALC が SALC 総数の 50% 以上を占める
- ホットパスが現状の `design_matrix_*` ではなく Monte Carlo sweep 風の
  「atom 単位で全 (l, m) を引く」アクセスパターンに変わる
- inner `sh_values` 確保の構造的撤廃（`MMatrix{N, 2L+1}` 化など）が許容される

それまでは Magesty 単独で **`MySphericalHarmonics` を維持**。SpheriCart は
SpinClusterMC 側でのみ採用、Magesty には依存追加不要。

### 残置物

- `test/component_test/test_sphericart_agreement.jl` は規約 drift 検出として
  **残存させる**（SpheriCart バージョン上げ時に Magesty 規約と乖離したら即検出）。
- `test/benchmark_optimize_sphericart.jl` / `test/benchmark_sphericart.jl` も保持。
  将来のアクセスパターン変更時に再評価できるベースラインとして有用。
- `Project.toml [compat]` の `SpheriCart = "0.2.3"` も `[extras]` 経由でテスト時のみ
  読み込まれる構成なのでそのまま。

---

## 設計案: Optimize.jl の estimator dispatch リファクタリング（2026-05-13、未着手）

**対象**: `src/Optimize.jl` の回帰手法ディスパッチと `elastic_net_regression`。

**動機**: 将来的に Lasso / Ridge / Bayesian / NNLS / 制約付き回帰などを追加する可能性がある。現状は型階層（`AbstractEstimator` / `OLS` / `ElasticNet`）は存在するが、`_fit_sce_model_internal` で `isa` 分岐しており、多重ディスパッチの利点を活かしていない。さらに `elastic_net_regression` が「重み付き問題組み立て」「求解」「j0 抽出」を 1 関数に混在させているため、新しい回帰を足すと組み立て・後処理がコピペで増殖する。

### 現状の構造的弱点

`elastic_net_regression`（L938–998）が独立した 3 責務を抱えている：

1. **問題組み立て**（estimator 非依存）— energy/torque を √weight でスケール、縦連結、bias 列整合。
2. **求解**（estimator 依存）— `X \ y` または `ridge(X, y, λ; bias=false)`。
3. **後処理**（estimator 非依存）— `j0 = mean(y_E - X_E[:,2:end] * jphi)` で bias 抽出。

`_fit_sce_model_internal` の `isa` 分岐（L837–857）は型を持ちながら dispatch しない「両方の悪いとこ取り」状態。

### 提案：3 層分離 + 多重ディスパッチ

```julia
abstract type AbstractEstimator end

struct OLS        <: AbstractEstimator end
struct Ridge      <: AbstractEstimator; lambda::Float64 end
struct Lasso      <: AbstractEstimator; lambda::Float64 end
struct ElasticNet <: AbstractEstimator; alpha::Float64; lambda::Float64 end
# struct NNLS, ConstrainedOLS, ElasticNetCV, ...

# 第1層: 問題組み立て（estimator 非依存）
assemble_weighted_problem(Xe, Xt, ye, yt, weight) -> (X, y, bias_col)

# 第2層: 求解（estimator ごとに dispatch、新しい回帰の追加点はここだけ）
solve_coefficients(::OLS,        X, y; bias_col) = X \ y
solve_coefficients(e::Ridge,     X, y; bias_col) = ridge_solve(X, y, e.lambda; bias_col)
solve_coefficients(e::Lasso,     X, y; bias_col) = lasso_solve(X, y, e.lambda; bias_col)
solve_coefficients(e::ElasticNet,X, y; bias_col) = enet_solve(X, y, e.alpha, e.lambda; bias_col)

# 第3層: 係数抽出（estimator 非依存）
extract_j0_jphi(j_values, Xe, ye) -> (j0, jphi)
```

これで `_fit_sce_model_internal` は 3 行に縮む：

```julia
X, y, bias_col = assemble_weighted_problem(Xe, Xt, ye, yt, weight)
j_values        = solve_coefficients(estimator, X, y; bias_col)
return extract_j0_jphi(j_values, Xe, ye)
```

**新しい回帰の追加 = 新しい struct + `solve_coefficients` の 1 メソッドだけ**（open/closed 原則）。`isa` 分岐の `else throw(ArgumentError(...))` も不要（未定義型は `MethodError` で自動失敗）。

### 設計上の判断ポイント

- **bias 列の扱い**：「bias を正則化から除外」は estimator 共通要件。`bias_col` を引数で渡すことで、各 solver が自分の方法で除外できる（自前実装なら `lambda_vec[bias_col]=0`、GLMNet.jl 系なら `intercept=true`）。今のうちにこのインターフェースを引き出しておく。
- **反復解法のオプション**（`max_iter`, `tol`, `warm_start`）: 各 estimator struct の field として持たせる。Holy traits は estimator が 10 を超えるまで導入しない。
- **CV / λ-path**：`ElasticNetCV` のような**別の estimator 型**として `solve_coefficients` を定義する。`ElasticNet` の中に CV ロジックを混ぜない。
- **外部パッケージ依存**（Lasso.jl, GLMNet.jl 等）：`solve_coefficients` の内部に閉じ込め、上位 API には漏らさない。

### 付随で整理したい点

- `ElasticNet.alpha` が "currently unused"（L46, L931）。L1 を入れる予定がないなら型名を `Ridge` に変えた方が誤解がない。L1 を入れる予定なら docstring に「現状 α≠0 でも無視される」を明記。
- `Optimizer` 構造体コンストラクタ（L94–104）で `alpha`, `lambda` を位置引数で取りつつ `estimator::AbstractEstimator = ElasticNet(alpha=alpha, lambda=lambda)` をデフォルトにしている設計は曖昧。estimator を渡したら位置引数が無視されるのか上書きされるのか呼び出し側から読み取れない。リファクタと合わせて `Optimizer(..., estimator, weight)` に寄せる。

### やらないこと（現時点）

- StatsAPI.jl / MLJ.jl 互換化（外部公開時に検討、今は内部 API で十分）。
- Holy traits（`CapabilityTrait` 系）— 4–5 種類なら多重 dispatch で捌ける。
- `AbstractEstimator{T}` のような抽象パラメータ化（YAGNI）。

### 進め方

API シグネチャに影響する中規模リファクタなので、着手前に `docs/specs/[YYMMDD]-estimator-dispatch/` を切って requirements.md / design.md / tasklist.md で合意する。実装単独で着手しない。

### 連動箇所

- `_fit_sce_model_internal`, `fit_sce_model_ols`, `fit_sce_model_elastic_net`, `elastic_net_regression`（L829–998）
- `Optimizer` 構造体の外側 constructor（L94–283）
- `fit_sce_model`（L678–790）の `estimator` 引数まわり
- `export` 一覧（L24）

---

## リファクタリング候補スイープ（2026-05-13、未着手）

**目的**: `src/` 全体（~9.4k 行）を対象に「構造的に整理すべき」候補を洗い出した結果。Optimize.jl は前セクションで個別カバー済みなので、本セクションでは触れない。

**着手判断**: 個別の項目を実装する前に CLAUDE.md の方針に従い、中規模リファクタリングは spec 化（`docs/specs/[YYMMDD]-[slug]/`）が必要。バグ修正やドキュメント追加で完結する小規模項目は spec なしで進めてよい。「要 verify」マーク付きの項目は、着手前に該当ファイルを読んで実態を確認する。

### 🔴 高優先度（API・数値結果・既知バグに関わるもの）

#### R1. `System` / `build_sce_basis` / `SpinCluster` の構築シーケンス重複

**対象**: `src/Magesty.jl`

`System(input_dict)`（L132–148）、`build_sce_basis`（L172–179）、`SpinCluster(input_dict)`（L286–307）の 3 箇所で、`Config4System → Structure → Symmetry → Cluster → BasisSet` の同じ 5 ステップが繰り返されている。さらに `build_sce_basis_from_xml`（L205–230）は BasisSet 部分だけ XML 読み込みに置換した 4 つ目のバリアント。

**改善案**: プライベートヘルパー `_build_core_components(config; verbosity, basisset_override=nothing)` を抽出して 4 箇所で再利用。`basisset_override` を渡せば XML 読み込みパスにも対応できる。

**連動箇所**: なし（モジュール内に閉じている）。

#### R2. `write_xml` の 3 シグネチャと位置/キーワード引数の混在

**対象**: `src/Magesty.jl` `write_xml`（L509–546 付近）

3 つのシグネチャ：
- `write_xml(sc::SpinCluster, filename, write_jphi=...)`
- `write_xml(system::System, filename)`
- `write_xml(structure, symmetry, basis_set, optimize, filename, write_jphi=...)`

`write_xml(system, "out.xml", false)` のような呼び出しで第 3 引数の意味が型ディスパッチによって変わり、ユーザーから読み取れない。

**改善案**: `write_xml(target; filename, write_jphi=true)` のキーワード引数中心 API に統一。

**注意**: 公開 API のシグネチャ変更なので、着手する場合は spec 化必須（CLAUDE.md「実装せず、必ず確認する」項目）。

#### R3. `BasisSets.jl:554` の TODO とコメントアウト処理（既知バグ）

**対象**: `src/BasisSets.jl` projection_matrix_coupled_basis 周辺

L554 に「Remove this exceptional handling when the bug is fixed in the projection matrix construction」のコメントとともに、コメントアウトされた `if sort(l_vec) == [1, 3]; continue; end` が残存。projection matrix 構築のバグ回避だったと思われるが、現在は完全にコメントアウト状態。

**改善案**: 以下のいずれか。
1. バグを特定して根本修正する（理想だが調査コストが読めない）。
2. 「known limitation」として再現条件・回避策・追跡 ID を本ファイルに正式に書き、ソースのコメントアウトコードは git 履歴に任せて削除する。

**連動箇所**: `test/component_test/` に該当ケース（`l_vec == [1, 3]`）のテストがあるか確認。

#### R4. `ConfigParser.jl` のデフォルト値・検証ロジックの 3 箇所散在

**対象**: `src/utils/ConfigParser.jl`

`DEFAULT_VALUES_SYSTEM` / `DEFAULT_VALUES_OPTIMIZE` Dict、`ValidationRule` struct、`Config4*.__init__` 内のハードコードされた `throw(ArgumentError(...))` チェックが分散している（要 verify）。新フィールドを追加すると 3 箇所編集が必要。

**改善案**: `ValidationRule` を inner constructor で declarative に apply、デフォルト値も同じテーブルから注入する。

**連動箇所**: `Config4System`, `Config4Optimize`、TOML スキーマドキュメント。

#### R5. `AngularMomentumCoupling.jl` の DP ロジック重複

**対象**: `src/utils/AngularMomentumCoupling.jl` `coeff_tensor_complex`（L137–212）と `coeff_tensor_complex_mindexed`（L279–359 付近、要 verify）

3 ステージの CG 漸化 DP（N ≥ 3）がほぼ同じアルゴリズムで 2 度実装されている。indexing 方式（`CartesianIndices` vs m-direct `OffsetArray`）のみ違う。

**改善案**: 共通の DP コアを 1 関数に抽出し、indexing 方式を wrapper で選択する。

**連動箇所**: `build_all_complex_bases`, `build_all_real_bases` が両方を呼んでいる。

### 🟡 中優先度（保守性・スケーラビリティ）

#### R6. `MySphericalHarmonics.jl` の buffered / non-buffered API の整理

**対象**: `src/utils/MySphericalHarmonics.jl`

`Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` で buffered 版と non-buffered 版が並列に存在。buffer サイズ要件（`length(buf) >= l - |m| + 1` 等）の docstring が散在しており、新規呼び出し時に必要サイズの根拠が辿りにくい。

**改善案**: buffer サイズ仕様を 1 箇所（例: モジュール冒頭の docstring または `BUFFER_SIZE_*` 定数）に集約し、各 buffered 関数の docstring から参照する。可能なら `@boundscheck` で buffer サイズの assertion を追加。

#### R7. `xml_io.jl` の XML タグ・属性のマジックストリング

**対象**: `src/utils/xml_io.jl`

`"Magesty"`, `"System"`, `"SALC"`, `"basis"` などのノード名・属性名、`"%.15e"` のフォーマット文字列が read/write 両側で hardcoded（CLAUDE.md「連動箇所: SCE 係数の入出力」で整合性が求められている部分）。

**改善案**: タグ・属性・フォーマット定数を 1 箇所に集約し、read/write が共有する。スキーマ変更が片側だけに入り込むリスクを排除。

**連動箇所**: BasisSet の XML ラウンドトリップテスト（あれば）。

#### R8. `SortedContainer.jl` の 3 コンテナの実装重複

**対象**: `src/common/SortedContainer.jl`

`SortedVector`, `SortedUniqueVector`, `SortedCountingUniqueVector` が `AbstractSortedVector` を継承する設計だが、`push!` / `append!` / `delete!` の本体ロジックが重複している（要 verify）。

**改善案**: 共通の sorted insertion / deletion 戦略を base に統合し、uniqueness / counting を直交した拡張として実装する。

#### R9. 公開 API の docstring 不揃い

**対象**: 主に `src/SpinConfigs.jl` `get_j0`, `get_jphi`, `get_j0_jphi`（L683–701 付近）と export 一覧全体

CLAUDE.md は「エクスポートされる API には明示的な型アノテーションと docstring を付ける」「Julia 標準形式（`# Arguments` / `# Returns` / `# Examples`）」を要求しているが、一部の getter で Arguments / Returns セクションが欠落。

**改善案**: `src/Magesty.jl` の export 一覧（L77–84）を起点に audit し、全 export 関数の docstring を Julia 標準形式に統一する。

#### R10. `Basis.jl` `CoupledBasis` constructor の不変条件検証の散在

**対象**: `src/types/Basis.jl:43–79`

docstring では `atoms` が "sorted" であるべきと書かれているが、inner constructor で assertion されていない（要 verify）。N=1 special case の handling、Lseq 長チェック、エラーメッセージのフォーマットが不揃い。

**改善案**: すべての不変条件を inner constructor 1 箇所に集約し、エラーメッセージのフォーマットを統一する。

#### R11. `EnergyTorque.jl` の `(4π)^(n_C/2)` スケーリングの重複と物理意味未注釈

**対象**: `src/utils/EnergyTorque.jl` `calc_energy`（L59 付近）, `calc_torque`（L134 付近）, および `src/Optimize.jl` `build_design_matrix_energy` / `build_design_matrix_torque` 内の同様の式

物理由来（technical_notes 参照）の docstring が欠落しているうえ、3 箇所以上に同じ式が散らばっている。スケーリングの規約が変わったら全箇所同期が必要だが、現状は静かに drift しうる。

**改善案**: 定数または `_sce_scaling_factor(n_C)` helper として 1 箇所に集約し、由来（[Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/)）を docstring に明記。

**連動箇所**: CLAUDE.md「連動箇所」セクションに追加する候補。

### 🟢 低優先度（参照のみ）

- **`BasisSets.jl` `construct_basislist`（L957–1000 付近）**: "backward compatibility" コメント付きで、内部呼び出しがあるかどうか要 verify。なければ削除。
- **`SpinConfigs.show`（L175–197）と `Clusters.print_cluster_stdout`（L506–566）**: 出力フォーマットの統一余地。`PrettyTables.jl` 検討の余地もあるが新規依存追加なので保留。
- **`AtomCells.jl` `AtomCell` と `AtomicIndices.jl` `SHSiteIndex`**: `isless` / `==` / `hash` の実装パターンが重複。**注意**: 物理意味が違うので無理に共通化するとバグの温床。慎重に。
- **`common/version.jl`**: hardcoded version 文字列。Project.toml から動的取得する仕組みの検討余地。
- **`Magesty.jl` `print_header`（L743–757）**: version / Julia info / threads / timestamp の 4 責務。verbose ロギングを整備する際に分割。

### Optimize.jl 関連項目（クロスリファレンス）

`alpha` パラメータの "currently unused"、`_fit_sce_model_internal` の `isa` 分岐、`elastic_net_regression` の責務混在、`Optimizer` コンストラクタの曖昧な API は、本セクションの**前にある「設計案: Optimize.jl の estimator dispatch リファクタリング」で詳細カバー済み**。重複記述しない。

---

## 未着手の改善案

### 🔴 候補F: `∂ᵢZlm_unsafe` での `P̄ₗₘ`/`dP̄ₗₘ` 統合計算（着手中）

**対象**: `src/utils/MySphericalHarmonics.jl` `∂ᵢZlm_unsafe(l, m, uvec, buf)`

**現状**: 候補E で alloc は解消したが、`∂ᵢZlm_unsafe` は `P̄ₗₘ(l, n, z, buf)` と
`dP̄ₗₘ_unsafe(l, n, z, buf)` を続けて呼ぶため、内部の `dnPl(x, l, n)` と
`dnPl(x, l, n+1)` がほぼ同じ Legendre キャッシュを 2 回構築している
（`_unsafednPl!` は毎回 `collectPl!` + `n` 段の漸化式を実行）。

**提案**: `LegendrePolynomials._unsafednPl!` を 1 回だけ呼んで am 階導関数の
キャッシュを作り、続けて漸化式を 1 ステップだけ手動で進めて (am+1) 階導関数も
取り出す内部ヘルパ `_legendre_pair_unsafe!(buf, x, l, am) -> (P_am_l, P_am1_l)` を追加。
`∂ᵢZlm_unsafe(l, m, uvec, buf)` 内でこれを呼び、`plm` と `dplm` を 1 回の
キャッシュ構築で得る。

**期待効果**:
- Legendre 部分が約 1.7〜2×
- `∂ᵢZlm_unsafe` 全体で **1.3〜1.5×**

**実装規模**: 中。`LegendrePolynomials` の内部関数 `_unsafednPl!` と `dPl_recursion`
（いずれも非エクスポート）を直接呼ぶため、上流の API 変更に対する追従が必要。
等価性テストで担保する。数値結果は不変。

**注意**:
- `_unsafednPl!` は cache を破壊的に書き換える。漸化式の 1 ステップ追加分は
  既存 `_unsafednPl!` の内部実装と整合させる
- buffer サイズ要件は変わらず `length(buf) >= l - |m| + 1`

### 🔴 候補E: `Zₗₘ_unsafe` のバッファ事前確保（Magesty.jl 側、最大の改善余地）

**対象**: `src/utils/MySphericalHarmonics.jl` `P̄ₗₘ`, `dP̄ₗₘ_unsafe`, `Zₗₘ_unsafe`, `∂ᵢZlm_unsafe`

**現状**: `_update_atom_zlm_cache!`（SpinClusterMC 側 sweep の 43%）の主要コストは
`LegendrePolynomials.dnPl(x, l, n)` が毎回 `zeros(Float64, l-n+1)` をヒープに確保していること。
1 sweep あたり 1024 回・29.7 KB のアロケーションが発生し、GC 圧の主因にもなっている。

**提案**: `LegendrePolynomials.dnPl` には事前確保バッファを受け取る `dnPl(x, l, n, A)` 版がある。
これを利用して以下の buffered API を追加する（既存 API は据え置き）:

- `P̄ₗₘ(l, m, r̂z, buf)` — `dnPl(x, l, |m|, buf)` を経由
- `dP̄ₗₘ_unsafe(l, m, r̂z, buf)` — `dnPl(x, l, |m|+1, buf)` を経由
- `Zₗₘ_unsafe(l, m, uvec, buf)` — 上記を経由
- `∂ᵢZlm_unsafe(l, m, uvec, buf)` — 上記を経由

呼び出し側（SpinClusterMC）は per-thread `Vector{Float64}(undef, max_l + 2)` を確保して渡す。

**期待効果**:
- `dnPl` 単体で 5.3×（16 ns → 3 ns）
- `_update_atom_zlm_cache!` で 2〜3× 高速化
- sweep 全体で **1.3〜1.5×**（57.9 μs → 約 40 μs）
- GC オーバーヘッド（5.1%）も同時に解消

**実装規模**: 中。数値結果は不変であるべき（等価性テストで担保）。

**注意**:
- buffer サイズ不足は `BoundsError`／無音破壊の原因になりうる。サイズ要件は `length(buf) >= l - |m| + 1`（dP̄ 側は `l - |m|`）。`max_l + 2` を確保し、`@boundscheck` で守るか先頭で 1 度 `@assert` する
- 並列化された呼び出し箇所では必ず per-thread buffer（共有すると race）
- 呼び出し側（SpinClusterMC）の修正は本リポジトリ側で別途実施

### 🟡 SH バッファの cache 引数化（旧 #1/#2）

**対象**: `src/Optimize.jl` `design_matrix_energy_element`, `calc_∇ₑu!`

`mutable struct SHCache` で `sh_values` / `atom_grad_values` をキャッシュ化し、`build_design_matrix_*` で per-thread に確保して引数で渡す案。一度試行したが、

- `calc_∇ₑu!` は `atom_site_idx == 0` で多くが早期 return するため、関数冒頭の eager allocation が逆効果
- SVector 版 (時間 -4.4% / メモリ +14%) と Vector 版 (時間 +6% / メモリ -3%) でトレードオフ対立
- BenchmarkTools 実測ノイズ ±5% 程度で有意性グレー

として見送り。再検討するなら:
- 早期 return パターンに合わせた lazy 確保（`Ref{Bool}` フラグ＋初回 work iter で確保）
- `cb1.atoms` を `SVector{N,Int}` 化して N を型パラメータにし、`MVector{N,Int}` でスタック確保
- `Vector{SVector{3,Float64}}` のメモリ計上が増える原因の特定（`@code_warntype` / `@allocated` で per-line 計測）

詳細経緯: `.claude/bench_log.md` の "#1/#2" セクション。

### 🟢 `cb1.atoms` の SVector 化（旧 #10 の本筋）

**対象**: `src/types/Basis.jl` `CoupledBasis`, `CoupledBasis_with_coefficient`

現在 `atoms::Vector{Int}`。`SVector{N,Int}` 化すれば `projection_matrix_coupled_basis` の `map(atom -> symmetry.map_sym[atom, n], cb1.atoms)` がゼロアロケーションになる。ただし `N` を型パラメータにする必要があり、Basis 構築側の API 変更も伴う。`#10` で簡易対応（ループバッファ hoist）したのでアロケーション影響は限定的だが、SVector 化すれば `find_translation_atoms` 等の引数型も整い、関連箇所の最適化余地が広がる。

### 🟢 候補（要検証、メモリ/allocs 主目的）

第二弾調査で試行→ B/F だけ採用。以下はメモリ/allocs では効果あったが時間は noise レベル（±5%）に収まっていたため非採用。GC 圧の高いシナリオなどで再検討余地あり。

- **`temp_projection_mat[row_range, col_range] = rot_mat * phase` の broadcast in-place 化**: BasisSet メモリ -3.9%、allocs -224k。1 行修正だが時間効果なし。
- **`tensor_inner_product` の `sum(conj.(t1) .* t2)` を融合 `@simd` ループに置換**: BasisSet メモリ -8.8%、allocs -224k。時間は noise。

詳細測定値は `.claude/bench_log.md` の "C" / "D" セクション参照（perf-2 ブランチ commit 履歴に残存）。

### スコープ外（やらない）

- 旧 #7 (BasisSets `round.()` / マスキング融合): ボトルネックでないため対応不要
- 旧 #8 (`eigenvec` スライス `@views` 化): ボトルネックでないため対応不要

### 教訓: 試したが逆効果だったアプローチ

第二弾実装で試行→ regression のため棄却。以下のパターンは Julia コンパイラが既に最適化しているケースが多く、安易な手書き置換はパフォーマンス低下を招く:

- **`coeff_tensor[idx_buf..., mf_idx]` の splat を strides 計算で線形インデックス化**: torque 時間 +34%、energy 時間 +1%。プロファイル overhead で hotspot に見えても、Julia の vararg dispatch は小規模 N (=2) で十分高速化されており、手書き strides の方が遅い。
- **`mf_grad_contribution .+= coeff_val * product_other .* grad_atom` の 3 成分手展開**: 時間 +38%、allocs +88%、メモリ +38%。broadcast `.+=` はコンパイラが in-place fuse するが、手書き `MVector[i] += x` は setindex! オーバーヘッドが大きく逆効果になる場合あり。

→ プロファイル sample count は時間に比例しない。改善前後で必ず `@time` (5+ 試行) で計測する。
