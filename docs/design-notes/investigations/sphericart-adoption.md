# SpheriCart.jl 採用の可否

**Status**: 完了（2026-05-11）

**結論**: Magesty.jl では SpheriCart を採用しない（SpinClusterMC とは逆の結論）。

## 背景

SpinClusterMC 側調査（同日付）で Magesty `Zₗₘ_unsafe` ≡ SpheriCart `Y_l^m_real`
が bit-exact 一致することは確認済み。Magesty 本体でも採用が有利か調査。

## 動作確認

- `make test-sphericart` → `SpheriCart vs TesseralHarmonics agreement` で
  700 / 700 PASS（`test/component_test/test_sphericart_agreement.jl`、lmax=4 全 (l,m)、7 単位ベクトル）
- 数値同一性は完全（SpinClusterMC 側調査と同じ）

## 性能実測（`make bench-optimize-sphericart` 系）

`test/benchmark_optimize_sphericart.jl` は v1（drop-in: per-atom `compute`）と
v2（per-spinconfig 事前計算キャッシュ）両方を実装済み。median 比較:

### fege_2x2x2（64 atoms, lmax = 1, 146 SALCs, 100 spinconfigs）

| 関数 | baseline (Magesty) | v1 (SpheriCart drop-in) | v2 (SpheriCart 事前計算) |
|---|---|---|---|
| `design_matrix_energy_element` | **54.0 μs** | 133 μs (2.46× slow) | 113 μs (2.09× slow) |
| `calc_∇ₑu` | **8.4 μs** | 23.5 μs (2.80× slow) | 21.5 μs (2.55× slow) |
| `build_design_matrix_torque`（1 config） | **208 ms** | — | 503 ms (2.42× slow) |

### fept_tetragonal_2x2x2（16 atoms, lmax = 2, 31 SALCs, 30 spinconfigs）

| 関数 | baseline | v1 | v2 |
|---|---|---|---|
| `design_matrix_energy_element` | **27.7 μs** | 2.40× slow | 2.07× slow |
| `calc_∇ₑu` | **7.7 μs** | 3.11× slow | 2.71× slow |
| `build_design_matrix_torque` | **6.56 ms** | — | 2.39× slow |

メモリも v2 は torque 行列で **785 MiB** vs baseline **293 MiB**（2.7×）。

## なぜ SpinClusterMC とは逆の結果か

アクセスパターンの違い:

- **SpinClusterMC**: sweep 中の atom cache は「ある atom について全 (l, m) up to L_max」を
  欲しい。SpheriCart の SVector 一括返しが access pattern と完全一致 → **10-20× 高速化**。
- **Magesty `design_matrix_energy_element` / `calc_∇ₑu`**: ループは
  `for (site_idx, atom) in translated_atoms; l = cbc.ls[site_idx]; for m_idx in 1:2l+1`
  と、**site ごとに固定の `l`** で 2l+1 個の m 値しか要らない。SpheriCart は
  L_max まで全 (l, m) 計算するため、その他の (l', m') 計算は完全に無駄。
  L_max=2 でも (L+1)² = 9 個のうち欲しいのは ~1-3 個程度 → SIMD 並列の節約より
  「不要 (l, m) を計算する分のオーバーヘッド」が勝る。

v2（事前計算キャッシュ）でも勝てない理由:
- per-spinconfig の `Y_cache, dY_cache` 構築は cheap だが、各 `itrans × site` 反復で
  `sh_values = [[Y_cache[...] for mi in 1:2l+1] for si in 1:N]` と **inner Vector を
  毎回ヒープ確保**する構造になっており、確保コストが Magesty 内蔵 `Zₗₘ_unsafe` の
  値計算コストを上回る。

## Magesty 側で SpheriCart を活かす条件（参考）

理論上 SpheriCart 有利になる構造は:
1. すべての SALC の `cbc.ls` 最大値 ≈ 個別 `l` 値の平均（i.e., 全 (l,m) を毎回使う）
2. inner loop の `Vector{Float64}` 確保を完全に排除（`MMatrix` / 線形 indexing）
3. spinconfig 単位の cache を `Matrix{Float64}` のまま row として読む

しかし現状 Magesty の SALC は cluster ごとに `ls` が混在（1〜lmax まで様々）、
高 lmax SALC は SALC 数全体の小さな割合。restructure コストに見合う改善は望めない。

## 採用条件下での将来再評価指針

以下が同時に成立する場合は再検証する価値あり:
- L_max ≥ 4 の SALC が SALC 総数の 50% 以上を占める
- ホットパスが現状の `design_matrix_*` ではなく Monte Carlo sweep 風の
  「atom 単位で全 (l, m) を引く」アクセスパターンに変わる
- inner `sh_values` 確保の構造的撤廃（`MMatrix{N, 2L+1}` 化など）が許容される

それまでは Magesty 単独で **`TesseralHarmonics` を維持**。SpheriCart は
SpinClusterMC 側でのみ採用、Magesty には依存追加不要。

## 残置物

- `test/component_test/test_sphericart_agreement.jl` は規約 drift 検出として
  **残存させる**（SpheriCart バージョン上げ時に Magesty 規約と乖離したら即検出）。
- `bench/benchmark_sphericart.jl` も保持。
  将来のアクセスパターン変更時に再評価できるベースラインとして有用。
- `Project.toml [compat]` の `SpheriCart = "0.2.3"` も `[extras]` 経由でテスト時のみ
  読み込まれる構成なのでそのまま。
