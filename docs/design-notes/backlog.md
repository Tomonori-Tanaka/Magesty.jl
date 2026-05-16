# パフォーマンス改善バックログ

軽量な思いつきメモ・第二弾調査の保留候補。実装に着手する場合は本ファイルから抽出して別 design-note または spec に昇格させる。

## ✅ 候補F: `∂ᵢZlm_unsafe` での `P̄ₗₘ`/`dP̄ₗₘ` 統合計算（完了）

**Status**: **完了** (commits `153e354` / `131d9a4`, 2026-05 頃)。
`src/TesseralHarmonics.jl` に private helper
`_legendre_pair_unsafe!(buf, x, l, am) -> (P_am_l, P_am1_l)` を導入し、
`_unsafednPl!` を 1 回だけ呼んで `am` 階のキャッシュを構築、続けて
`dPl_recursion` で 1 ステップだけ in-place に進めて `(am+1)` 階を取り出す形に。
buffered `∂ᵢZlm_unsafe(l, m, uvec, buf)` がこの helper を経由するようになり、
Legendre キャッシュ構築が 2 回 → 1 回に統合された。

**実装上の注意点（残し）**:
- `LegendrePolynomials._unsafednPl!` と `dPl_recursion` はいずれも非エクスポート
  内部関数のため、上流の API 変更に対する追従が必要。
- buffer サイズ要件はモジュール docstring の "Buffer requirements" に集約済み
  (`length(buf) >= l - |m| + 1`)。`@boundscheck checkbounds(buf, _required_buf_size(l, m))`
  で守られる。
- 数値結果は不変（等価性テストで担保）。

## ✅ 候補E: `Zₗₘ_unsafe` のバッファ事前確保（Magesty.jl 側完了）

**Status**: **Magesty.jl 側完了** (commit `8a4a17d`, 2026-05-11)。
`P̄ₗₘ` / `dP̄ₗₘ_unsafe` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` に 4 引数の buffered
オーバーロードを追加し、`LegendrePolynomials.dnPl(x, l, n, A)` を経由して
per-call `zeros(l-|m|+1)` の heap 確保を解消。既存 3 引数 API は据え置き。
buffer サイズ要件は `length(buf) >= l - |m| + 1`（モジュール docstring に集約、
候補F 完了時の `131d9a4` で `@boundscheck checkbounds(buf, _required_buf_size(l, m))`
にて担保）。

**ベンチ**: `Zₗₘ_unsafe` (l=4, m=2) 単体で 35.2 ns → 20.8 ns (1.70×, 80 B → 0)。
数値結果は等価（l=0..8 × 全 m × 9 方向の等価性テスト + `@allocated == 0`）。

**残作業（本リポジトリ外）**: caller 側 (`SpinClusterMC` / `JPhiMagestyCarlo.jl`)
で per-thread `Vector{Float64}(undef, max_l + 1)` を確保して buffered overload
を呼ぶ修正は別 repo で実施する（並列化箇所では per-thread buffer 必須）。

## 🟡 SH バッファの cache 引数化（旧 #1/#2）

**対象**: `src/Fitting.jl` `design_matrix_energy_element`, `calc_∇ₑu!`

`mutable struct SHCache` で `sh_values` / `atom_grad_values` をキャッシュ化し、`build_design_matrix_*` で per-thread に確保して引数で渡す案。一度試行したが:

- `calc_∇ₑu!` は `atom_site_idx == 0` で多くが早期 return するため、関数冒頭の eager allocation が逆効果
- SVector 版 (時間 -4.4% / メモリ +14%) と Vector 版 (時間 +6% / メモリ -3%) でトレードオフ対立
- BenchmarkTools 実測ノイズ ±5% 程度で有意性グレー

として見送り。再検討するなら:
- 早期 return パターンに合わせた lazy 確保（`Ref{Bool}` フラグ＋初回 work iter で確保）
- `cb1.atoms` を `SVector{N,Int}` 化して N を型パラメータにし、`MVector{N,Int}` でスタック確保
- `Vector{SVector{3,Float64}}` のメモリ計上が増える原因の特定（`@code_warntype` / `@allocated` で per-line 計測）

詳細経緯: `.claude/bench_log.md` の "#1/#2" セクション。

## 🔴 `cb1.atoms` の SVector 化（旧 #10 の本筋）— ruled out

**Status**: 2026-05-16 に試行し abandon。spec: `docs/specs/260516-coupled-basis-atoms-svector/`。

`CoupledBasis{R, N}.atoms` を `SVector{N, Int}` 化して再ベンチ。
`projection_matrix_coupled_basis` を含む SALC build が **+12% time / +6 MB / +140K allocs** 悪化、design-matrix 経路は変化なし。

**原因**: `CoupledBasis` は `SortedCounter{CoupledBasis}` / `Vector{CoupledBasis}` という UnionAll コンテナで保持される。イテレーション時に `cb.atoms` の型パラメータ `N` が型消去されるため、`SVector{N, Int}` のインデクシングが specialize されず、`Vector{Int}` (element 型固定・長さ動的) より遅くなる。design-matrix 側はホットパスが `where {R, N}` の function barrier で specialize されるので影響を受けない。

**再着手の前提条件**:

- コンテナの concrete 化（(R, N) ごとに分離した `Vector{CoupledBasis{R, N}}`、または function barrier による grouped iteration）
- `map(closure, atoms)` の closure capture を `let` 束縛で回避

これらが揃わない限り net negative。spec の "Outcome" セクション参照。

## 🟢 候補（要検証、メモリ/allocs 主目的）

第二弾調査で試行→ B/F だけ採用。以下はメモリ/allocs では効果あったが時間は noise レベル（±5%）に収まっていたため非採用。GC 圧の高いシナリオなどで再検討余地あり。

- **`representation_mat[row_range, col_range] = rot_mat * phase` の broadcast in-place 化**: SALCBasis メモリ -3.9%、allocs -224k。1 行修正だが時間効果なし。（旧名 `temp_projection_mat`、2026-05-16 改名）
- ~~**`tensor_inner_product` の `sum(conj.(t1) .* t2)` を融合 `@simd` ループに置換**: SALCBasis メモリ -8.8%、allocs -224k。時間は noise。~~ *(2026-05-16: 決定論的なスカラーループに置換済み。動機は cross-platform 一致 (BLAS/LAPACK SIMD lane width に依らない reduction order) で、メモリ/allocs 削減は副産物。`@simd` ではなくプレーンな `@inbounds for` を採用。)*

詳細測定値は `.claude/bench_log.md` の "C" / "D" セクション参照（perf-2 ブランチ commit 履歴に残存）。

## スコープ外（やらない）

- 旧 #7 (SALCBases `round.()` / マスキング融合): ボトルネックでないため対応不要
- 旧 #8 (`eigenvec` スライス `@views` 化): ボトルネックでないため対応不要

## 教訓: 試したが逆効果だったアプローチ

第二弾実装で試行→ regression のため棄却。以下のパターンは Julia コンパイラが既に最適化しているケースが多く、安易な手書き置換はパフォーマンス低下を招く:

- **`coeff_tensor[idx_buf..., mf_idx]` の splat を strides 計算で線形インデックス化**: torque 時間 +34%、energy 時間 +1%。プロファイル overhead で hotspot に見えても、Julia の vararg dispatch は小規模 N (=2) で十分高速化されており、手書き strides の方が遅い。
- **`mf_grad_contribution .+= coeff_val * product_other .* grad_atom` の 3 成分手展開**: 時間 +38%、allocs +88%、メモリ +38%。broadcast `.+=` はコンパイラが in-place fuse するが、手書き `MVector[i] += x` は setindex! オーバーヘッドが大きく逆効果になる場合あり。

→ プロファイル sample count は時間に比例しない。改善前後で必ず `@time` (5+ 試行) で計測する。
