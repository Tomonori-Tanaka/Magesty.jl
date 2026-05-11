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
