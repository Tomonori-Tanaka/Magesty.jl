# Bench Log

`design_note.md` の改善案を実装する際の前後ベンチマーク記録。

## 計測環境

- マシン: Darwin 24.6.0 (arm64)
- Julia: 1.12.6
- スレッド数: 4 (`--threads=auto`)
- ターゲット: `test/examples/fege_2x2x2/input.toml`
- スクリプト: `test/benchmark_basisset_hotspots.jl --samples 20 --evals 1 --profile-iters-* 0`

## ベースライン (commit `0acafc9` 直後 / perf branch HEAD)

| Function | Time median | Memory | Allocs |
|---|---|---|---|
| `listup_coupled_basislist` | 1.56 μs | 2.08 KiB | 49 |
| `projection_matrix_coupled_basis` | **79.5 ms** | **208.7 MiB** | **1,978,133** |
| `BasisSet` constructor | 149.8 ms | 730.9 MiB | 11,890,219 |

---

## #4: SALC 投影行列のバッファ再利用

**修正対象**: `src/BasisSets.jl` `projection_matrix_coupled_basis`

**経緯**: 当初は `@threads :static` 並列化を試みたが、外側 `_compute_salc_groups` の `@threads`（BasisSets.jl:179）からも呼ばれており Julia の「`@threads :static` cannot be used concurrently or nested」制約に抵触。外側で既にスピン構成単位の並列化が効いているため、内側はシリアル維持しつつ `temp_projection_mat` をループ外で 1 回確保 → `fill!` で再利用する形に変更。

### Before / After (`projection_matrix_coupled_basis`)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 79.5 ms | **60.4 ms** | **-24 %** |
| Memory | 208.7 MiB | **77.0 MiB** | **-63 %** |
| Allocs | 1,978,133 | 1,976,792 | -0.07 % |

### Before / After (`BasisSet` constructor)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 149.8 ms | 149.6 ms | -0.1 % |
| Memory | 730.9 MiB | **491.4 MiB** | **-33 %** |
| Allocs | 11,890,219 | 11,825,848 | -0.5 % |

**所感**: 内側関数で 24% の高速化＋メモリ -63%。BasisSet 全体での時間短縮はほぼ無いが、メモリ消費 -33% は GC 圧低減として効く（並列化時にスケールしやすくなる）。残る allocations は内側ループの `atoms_shifted_list` comprehension などが支配的（→ 後続の #10 で対処予定）。

---

## #9: Wigner D 行列の precompute

**修正対象**: `src/BasisSets.jl` `projection_matrix_coupled_basis`

`base_rot_mat = Δl(Lf, α, β, γ)` は `n` (symop) のみに依存し、`time_rev_sym` には依存しない。元コードでは `for ..., time_rev_sym in [false, true]` の二重展開により 2×nsym 回計算していた。`Vector{Matrix{Float64}}` (length nsym) として関数冒頭で precompute し、ループ内では参照のみ。

### Before / After (`projection_matrix_coupled_basis`)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 60.4 ms | 62.3 ms | +3 % (noise) |
| Memory | 77.0 MiB | 76.7 MiB | -0.4 % |
| Allocs | 1,976,792 | 1,968,762 | -0.4 % |

### Before / After (`BasisSet` constructor)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 149.6 ms | **90.2 ms** | **-40 %** |
| Memory | 491.4 MiB | **380.7 MiB** | **-23 %** |
| Allocs | 11,825,848 | **8,757,988** | **-26 %** |

**所感**: 単発の `projection_matrix_coupled_basis` 計測ではほぼ差が出ないが、`_compute_salc_groups` が key group 単位で多数呼ばれ、その内部で Δl の重複計算が積み上がっていたため、BasisSet 全体では -40% 時間短縮と allocations -3M 規模の改善。`Δl` 内部で行列確保しているのが主因と推察。

---

## #6: `build_design_matrix_torque` の preallocation 化

**修正対象**: `src/Optimize.jl` `build_design_matrix_torque`

`design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)` に各スレッドのブロック行列を入れて最後に `vcat(...)` していた。これを `Matrix{Float64}(undef, num_spinconfigs * 3 * num_atoms, num_salcs)` で 1 度確保し、各スレッドが `sc_idx` の disjoint な行レンジに直接書き込む形に変更。

### Before / After (`build_design_matrix_torque`, fege 20 spinconfigs)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 1.798 s | 1.822 s | +1 % (noise) |
| Memory | 6.46 GiB | 6.45 GiB | -0.15 % |
| Allocs | 240,384,112 | 240,384,047 | -65 allocs |

**所感**: 単独効果は微小。総アロケーション量は `calc_∇ₑu!` の内部割り当て (#2) が支配的なため、#6 の vcat 廃止だけでは数値が動かない。ただし以下の点で価値あり:
- 中間 `Vector{Matrix}` と最終 `vcat` の中間コピーを廃止（構造的にクリーン）
- #2 (calc_∇ₑu! のバッファ前確保) を入れた後に、相対的な寄与が見える可能性

---

## #5: `construct_map_sym` の SVector 化

**修正対象**: `src/Symmetries.jl` `construct_map_sym`

`@threads` ループ内の `rotation*x_frac[:, iat]` / `x_frac[:, jat] - local_x_new` / `abs.(diff) .% 1.0` を SVector ベースに置換。`MVector{3}` スクラッチ (`local_x_new`, `local_tmp`) を撤去し、回転行列・並進ベクトルも `SMatrix` / `SVector` で 1 度確保。距離判定を二乗距離ベースにして `norm` の `sqrt` も省略。

### Before / After (`Magesty.System(input)`, fege)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 152.4 ms | 153.6 ms | +0.8 % (noise) |
| Memory | 716.1 MiB | 708.1 MiB | -1.1 % |
| Allocs | 17,162,253 | 16,953,261 | **-209k** |

**所感**: `construct_map_sym` は System 構築 1 回限りの呼び出しで、System 全体に占める割合は小さいため目立たないが、内側ループのヒープ割り当てを完全に消した（SVector がスタックに乗る）。allocs -209k はそのまま GC 負荷低減として効く。

---

## #11: `design_matrix_energy_element` の `@inbounds` 補完

**修正対象**: `src/Optimize.jl` `design_matrix_energy_element`

`calc_∇ₑu!` の `for itrans` には既に `@inbounds` が付いていたが、エネルギー側 `design_matrix_energy_element` の `for itrans in symmetry.symnum_translation` (line 355) に欠けていた。これを `@inbounds for itrans` に変更し、ループ全体（map_sym indexing、`cbc.ls[]`、`sh_values[][]`、`coeff_tensor[idx_buf..., mf_idx]` など）の境界チェックを抑制。

### Before / After (`build_design_matrix_energy`, fege 20 spinconfigs)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 1.483 s | **1.458 s** | **-1.7 %** |
| Memory | 1.81 GiB | 1.81 GiB | 0 % |
| Allocs | 91,976,323 | 91,976,323 | 0 |

**所感**: 境界チェック除去のみなので allocs は不変。1.7% の高速化はそのまま CPU 命令削減効果。残りの allocations 削減は #1 が必要。

---

## #3: `idx_buf` / `other_sites` / `other_dims` の hoist + comprehension 撤去

**修正対象**: `src/Optimize.jl` `design_matrix_energy_element`, `calc_∇ₑu!`

両関数の `for itrans` ループ内で毎回再構築されていた以下のバッファをループ外へ hoist:
- `idx_buf::Vector{Int}` （N 要素、書き換えのみ）
- `other_sites_buf::Vector{Int}` / `other_dims_buf::Vector{Int}` （N-1 要素）

特に `calc_∇ₑu!` の `other_sites = [s for s in 1:N if s != atom_site_idx]` は comprehension で毎回 Vector 確保していたが、`atom_site_idx` 単純比較ループに置換し、事前確保バッファに書き込む形へ。`design_matrix_energy_element` 側は `atom_site_idx` 依存ではないので関数冒頭で 1 度確保するだけ。

### Before / After (`build_design_matrix_energy`, fege 20 spinconfigs)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 1.529 s | **1.436 s** | **-6.1 %** |
| Memory | 1.81 GiB | 1.77 GiB | -2.2 % |
| Allocs | 91,976,323 | 90,631,363 | **-1.34M** |

### Before / After (`build_design_matrix_torque`, fege 20 spinconfigs)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 1.765 s | **1.262 s** | **-28.5 %** |
| Memory | 6.45 GiB | 5.73 GiB | **-11.2 %** |
| Allocs | 240,384,047 | 206,804,527 | **-33.6M** |

**所感**: torque 側の comprehension が支配的なアロケーション源だったため大幅改善。energy 側は元から `1:(N-1)` の UnitRange を使っており影響は限定的だが、`dims[other_sites]` Vector 確保と `idx_buf` の per-iter 確保が消えて 1.3M allocs 削減。fege の integration test 全体時間が 65s → 52s に短縮。
