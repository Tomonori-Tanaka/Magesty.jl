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

---

## #1/#2: SH バッファの cache 引数化（**スキップ**）

**修正対象（試行）**: `src/Optimize.jl` `design_matrix_energy_element`, `calc_∇ₑu!`, `build_design_matrix_energy`, `build_design_matrix_torque`

design_note 通り `mutable struct SHCache` を導入し、`build_design_matrix_*` で per-thread キャッシュを構築 → calc 関数に渡す形を試した。

### 試行内容
- `SHCache` 型（`sh_values::Vector{Vector{Float64}}` + `atom_grad_values`）を定義
- `ensure_capacity!` で逐次 grow
- 関数シグネチャに `cache::SHCache = SHCache()` を追加（位置引数）
- `build_design_matrix_torque` で `caches = [SHCache() for _ in 1:Threads.maxthreadid()]` を確保し、`Threads.threadid()` でインデックス

### 結果（fege 20 spinconfigs、3 構成比較）

| 構成 | Time | Memory | Allocs |
|---|---|---|---|
| キャッシュ無し (#3 baseline) | 1.401 s | 5.73 GiB | 206.8M |
| SVector キャッシュ (`Vector{SVector{3,Float64}}`) | 1.340 s (**-4.4%**) | 6.53 GiB (**+14%**) | 217.7M (+5.3%) |
| Vector キャッシュ (`Vector{Vector{Float64}}`) | 1.491 s (+6.4%) | 5.57 GiB (-3%) | 202.7M (-2%) |

### スキップ判断
- **トレードオフが対立**: 時間 ↔ メモリでどちらか一方を犠牲にする形になり、決定的に勝つ構成がない
- BenchmarkTools 実測ノイズも ±5% 程度あり、時間差が有意かグレーゾーン
- API 変更（cache 引数追加）のコストに見合う改善が得られない
- `calc_∇ₑu!` は `atom_site_idx == 0` で多くが早期 return するため、cache を渡してもバッファ確保の節約効果が薄い
- 単純な hoisting 試行は早期 return パターンと相性が悪く悪化（事前に検証済み）

### 教訓
- 早期 return が支配的な関数では、関数冒頭の eager allocation は害になる
- mutable struct field アクセスと local Vector アクセスでメモリ計上が異なる挙動を示すケースがあり、ベンチで素直な結果が出ない
- `Vector{SVector{N,T}}` は要素が bitstype でも、確保サイズが OLD `Vector{Vector{Float64}}` より大きくなり全体メモリを押し上げる場合がある

#3 (idx_buf/other_sites の hoist) で既に大半のアロケーションは削減済み（torque -33.6M）なので、ここで打ち切って次の項目へ。

---

## #10: `atoms_shifted_list` comprehension 撤去

**修正対象**: `src/BasisSets.jl` `projection_matrix_coupled_basis`

ループ内側 `atoms_shifted_list = [symmetry.map_sym[atom, n] for atom in cb1.atoms]` の comprehension が (n, time_rev) × cb1 の組合せ分（fege で約 200 × 192 ≈ 38k）走り、毎回 Vector を確保していた。`coupled_basislist` 内の全 cb は同じクラスタサイズなので、関数冒頭で `atoms_shifted_list = Vector{Int}(undef, N_atoms)` を 1 度確保し、`@inbounds` で in-place 書き込みに変更。

### Before / After (`projection_matrix_coupled_basis`)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 61.6 ms | 62.0 ms | ±0 (noise) |
| Memory | 76.74 MiB | 73.75 MiB | -3.9 % |
| Allocs | 1,968,762 | 1,890,428 | **-78k** |

### Before / After (`BasisSet` constructor)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 89.9 ms | 92.0 ms | +2 % (noise) |
| Memory | 380.71 MiB | 372.16 MiB | -2.2 % |
| Allocs | 8,757,988 | 8,533,840 | **-224k** |

**所感**: 時間影響はないが、内側ループの隠れアロケーションを除去。コードもシンプル化（comprehension → 明示的ループ + `@inbounds`）。`cb1.atoms` を `SVector{N,Int}` 化すれば map で zero-alloc にできるが Basis 型の構造変更が必要なので深追いせず。

---

# 第二弾改善

第一弾マージ後に再度プロファイラで現状調査。残るホットスポットは:
- `coeff_tensor[idx_buf..., mf_idx]` の splat (energy line 410, torque line 620)
- `mf_grad_contribution .+= coeff_val * product_other .* grad_atom` (torque line 621)
- `build_design_matrix_energy` が単スレッド (utilization 50%)

## A: coeff_tensor の線形インデックス化（**スキップ**）

**修正対象（試行）**: `src/Optimize.jl` `design_matrix_energy_element`, `calc_∇ₑu!`

`idx_buf...` splat を strides 計算ベースの線形インデックスに置換。Vector 静的型でない splat は実行時 vararg を生成しヒープに確保するためプロファイル上は最大ホットスポットだった。

### 結果（fege 20 spinconfigs、10 trials @time）

| 関数 | Before | After | Δ |
|---|---|---|---|
| `build_design_matrix_energy` 時間 | 2.13 s | 2.16 s | +1 % |
| `build_design_matrix_energy` allocs | 90.6M | 72.5M | **-20 %** |
| `build_design_matrix_energy` メモリ | 1.77 GiB | 1.23 GiB | **-30 %** |
| `build_design_matrix_torque` 時間 | 1.36 s | 1.82 s | **+34 %** |
| `build_design_matrix_torque` allocs | 206.8M | 170.6M | -17.5 % |
| `build_design_matrix_torque` メモリ | 5.73 GiB | 4.65 GiB | -19 % |

**スキップ判断**: torque で時間 +34% の明確な regression。allocs/メモリは改善するが、Julia の splat dispatch が小規模 N (=2) で十分高速化されており、手書き strides の方が遅い結果に。コンパイラ最適化の余地は探ったが、`@inbounds` 追加・`strides_buf[last_site]` の hoist 等いずれも効果なし。リバート。

教訓: プロファイラの sample count は時間ではなくサンプリング頻度であり、低レベル splat はサンプル多くてもトータルでは早い場合がある。

## B: `build_design_matrix_energy` の `@threads` 並列化

**修正対象**: `src/Optimize.jl:300`

プロファイル utilization 50% (4 スレッド環境で実質 1 スレッド) を解消。`for i = 1:num_salcs` を `@threads` 化、各スレッドは disjoint な列 `design_matrix[:, i+1]` に書き込み。

### Before / After (`build_design_matrix_energy`, fege 20 spinconfigs、10 trials @time)

| Metric | Before | After | Δ |
|---|---|---|---|
| Time median | 2.13 s | **0.54 s** | **-75 %** (約 4x speedup) |
| Memory | 1.77 GiB | 1.77 GiB | 0 |
| Allocs | 90.63M | 90.63M | 0 |

FeGe integration test: **1m20.8s → 56.9s (-30%)**。

**所感**: 単純な並列化だが効果絶大。energy 側は既に hot spot は計算量で支配されていたため、4 スレッドでほぼ理想的にスケール。
