# Design Note

改善案・設計メモを記録するファイル。実装済みの履歴は `.claude/bench_log.md` と `git log` を参照。

---

## 未着手の改善案

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

### スコープ外（やらない）

- 旧 #7 (BasisSets `round.()` / マスキング融合): ボトルネックでないため対応不要
- 旧 #8 (`eigenvec` スライス `@views` 化): ボトルネックでないため対応不要
