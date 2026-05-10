# Design Note

改善案・設計メモを記録するファイル。

---

## パフォーマンス改善案（2026-05-10 調査）

`src/Optimize.jl` / `src/BasisSets.jl` / `src/Symmetries.jl` のホットパスを精査した結果。
すべての項目はスレッド並列との両立性を検討済み。

### 🔴 HIGH

1. **`Optimize.jl:371-380` `design_matrix_energy_element` の `sh_values` 再確保**
   - 並進ループ内で `sh_values = Vector{Vector{Float64}}(undef, N)` とその要素 `Vector{Float64}` を毎回確保。
   - 修正案: `cache::SHCache`（`sh_values_buf` フラット化版 + サイトオフセット表）を引数化。`build_design_matrix_energy` の `@threads` ループ直前で `[SHCache() for _ in 1:nthreads()]` を確保し、`cache_arr[threadid()]` で渡す。
   - 並列両立性: ✅ per-thread バッファで完全両立。

2. **`Optimize.jl:553-573` `calc_∇ₑu!` の `sh_values` / `atom_grad_values` 再確保**
   - 同型の問題。`build_design_matrix_torque` から `num_spinconfigs × num_atoms × num_salcs × num_cbcs × num_translations` 回呼ばれる。
   - 修正案: `cache::SHCache`（`sh_values_buf`, `grad_values_buf`, `idx_buf`）を引数化。呼び出し側 `@threads for sc_idx` 直前に per-thread で確保。
   - 並列両立性: ✅ per-thread バッファで完全両立。

3. **`Optimize.jl:579, 390, 582` `idx_buf` / `other_sites` / `CartesianIndices` の再構築**
   - `other_sites = [s for s in 1:N if s != atom_site_idx]` は新規 Vector。
   - 修正案: クラスタサイズ N が小さいため `NTuple{N-1,Int}` / `MVector{N,Int}` 化してスタック割り当てへ。
   - 並列両立性: ✅ スタック確保なので各タスク独立、共有状態なし。

4. **`BasisSets.jl:744-779` SALC 投影行列ループの `@threads` 並列化 + per-thread バッファ**
   - 現状はシリアル。`2*nsym` 回ごとに `temp_projection_mat = zeros(...)` を確保。
   - 修正案:
     - `(n, symop, time_rev_sym)` を平坦化したリストを `@threads` で分散。
     - スレッド数ぶんの `local_temp_mats` / `local_proj_mats` を確保し `fill!` で再利用。各スレッドはローカル `proj_mat` に累積し、最後に reduction sum。
     - `threadid()` ベースは Julia 1.9+ で不安定なので `@threads :static` 固定 or `ChunkSplitters.jl` / `OhMyThreads.jl` で分割。
   - 期待効果: 並列化で `nthreads` 倍程度の高速化 + per-thread バッファ再利用で確保コストもほぼゼロに。
   - 並列両立性: ✅ per-thread accumulator + reduction で完全両立。

5. **`Symmetries.jl:268-274` `diff = supercell.x_frac[:, jat] - local_x_new`**
   - 既存 `@threads` ループ内で 3 要素 Vector を毎回確保。
   - 修正案: `SVector{3,Float64}` 化（`@SVector` リテラル or 3 要素ループ展開）。
   - 並列両立性: ✅ スタック確保。

### 🟡 MEDIUM

6. **`Optimize.jl:487` `vcat(design_matrix_list...)` の大規模割り当て**
   - 結果行列 `(3*num_atoms*num_spinconfigs) × num_salcs` を最初に確保し、各スレッドが `sc_idx` の disjoint な行ブロックへ `@view` で書き込む。
   - 並列両立性: ✅ 行ブロックが disjoint なので競合なし。

7. **`BasisSets.jl:87-88, 100-101` `round.()` / マスキングのブロードキャスト連鎖**
   - 中間配列が複数生成。
   - 修正案: 融合ループ `for i in eachindex(eigenvecs); v=eigenvecs[i]; eigenvecs[i] = abs(v) < 1e-8 ? 0.0 : round(v, digits=10); end`。
   - 並列両立性: ✅ `_compute_salc_groups` の local 配列に対する操作。

8. **`BasisSets.jl:108-110` `eigenvec[coeff_start:coeff_end]` のスライスコピー**
   - ゼロチェック前は `@views` で参照し、push 時のみ `copy`。`norm(coefficient)` は `maximum(abs, view)` に置換。
   - 並列両立性: ✅ local データへの read。

9. **`BasisSets.jl:750-751` Wigner D 行列 `Δl(Lf, α, β, γ)` の precompute**
   - `(Lf, n, time_rev)` で完全に決まるため、`projection_matrix_coupled_basis` 冒頭（並列ループ前）で `Vector{Matrix{Float64}}`（`(n, time_rev)` でインデックス）を precompute し、ループ内は read-only。
   - ⚠️ `Dict` の lazy キャッシュ方式はスレッド書き込み競合があるため採用しない。事前計算限定。
   - 並列両立性: ✅ read-only 共有データなので安全。

### 🟢 LOW

10. **`BasisSets.jl:755` `atoms_shifted_list` の comprehension**
    - `cb1.atoms` が `SVector` 化されていれば `map(...)` で `SVector` 結果を返しゼロアロケーション化。`Structures.jl` 側との整合確認要。
    - 並列両立性: ✅ スタック確保。

11. **`Optimize.jl` 内側ループの `@inbounds` 漏れ**
    - 378 / 408 / 566 などの多重 indexing 行を `@inbounds begin ... end` で覆う。
    - 並列両立性: ✅ スレッディングと無関係。

### 推奨実装順

1. #4 SALC 投影行列の `@threads` 並列化 + per-thread バッファ
2. #9 Wigner D 行列の precompute
3. #1 / #2 / #3 `Optimize.jl` バッファ前確保（`cache` 引数追加リファクタを伴う）
4. #6 design_matrix の preallocation
5. #5 / #7 / #8（局所的最適化）
6. #10 / #11（クリーンアップ）

各ステップで `make test-all` / `make test-jet` を実行して回帰確認。代表ケースで `@btime` 比較し allocation 削減を計測する。
