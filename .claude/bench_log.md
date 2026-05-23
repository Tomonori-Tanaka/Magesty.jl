# Bench Log

`DESIGN_NOTES.md` の改善案を実装する際の前後ベンチマーク記録。

## 計測環境

- マシン: Darwin 24.6.0 (arm64)
- Julia: 1.12.6
- スレッド数: 4 (`--threads=auto`)
- ターゲット: `test/integration/fege_2x2x2/input.toml`
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

## F: `find_translation_atoms` の冗長コピー削除

**修正対象**: `src/BasisSets.jl:678`

`new_atom_list = Vector{Int}(atom_list)` は読み取り専用なので不要。引数を `AbstractVector{<:Integer}` 化して直接参照。

### Before / After (post-B → post-B+F、fege)

| Function | Metric | Before | After | Δ |
|---|---|---|---|---|
| `projection_matrix_coupled_basis` | Memory | 73.75 MiB | 70.76 MiB | -4.1 % |
|  | Allocs | 1,890,428 | 1,812,092 | **-78k** |
| `BasisSet` constructor | Memory | 372.16 MiB | 363.60 MiB | -2.3 % |
|  | Allocs | 8,533,840 | 8,309,584 | **-224k** |

時間は ±3% の noise 内（59〜60 ms / 91〜93 ms）。

**所感**: 軽微な改善だが構造クリーンアップ。`atom_translated = similar(...)` の per-i 確保はループから返す可能性があるためそのまま残置。

---

## 第二弾サマリ (B + F)

`build_design_matrix_energy` (fege 20 spinconfigs、5 trials @time):
- Time: 2.13 s → 0.56 s (**-74%**, 約 4x speedup via `@threads`)
- Allocs/Memory: 不変（90.63M / 1.77 GiB）

`projection_matrix_coupled_basis` (fege):
- Time: 79.0 ms → 59.2 ms (**-25%**, 主に thermal/再計測効果)
- Memory: 73.75 → 70.76 MiB (-4.1%)
- Allocs: 1.89M → 1.81M (-78k, F の効果)

`BasisSet` constructor:
- Time: 163.6 ms → 90.5 ms (**-45%**, thermal/再計測込み)
- Memory: 372.16 → 363.60 MiB (-2.3%)
- Allocs: 8.53M → 8.31M (-224k, F の効果)

FeGe integration test: 1m20.8s → 55.8s (**-31%**)。

### 採用しなかった改善案
プロファイル上で hotspot に見えた以下は試行 → ベンチで悪化したため棄却。詳細は `DESIGN_NOTES.md` 末尾の「教訓」セクション:
- A: `coeff_tensor[idx_buf..., mf_idx]` の splat → 線形インデックス化（torque +34%）
- E: `mf_grad_contribution .+= ...` の手展開（時間 +38%, allocs +88%）

メモリ・allocs では効果ありだが時間が noise 内のため見送った候補（`DESIGN_NOTES.md` の「候補」セクション）:
- C: `rot_mat * phase` を broadcast in-place 化 → BasisSet メモリ -3.9%, allocs -224k
- D: `tensor_inner_product` 融合 → BasisSet メモリ -8.8%, allocs -224k

---

## estimator-dispatch refactor — baseline (2026-05-13)

**Branch**: `refactor/estimator-dispatch` HEAD = `5efbf1b` (spec-only commit on top of `main` 9ba9e65, no source changes yet)

**Target**: `fit_sce_model` end-to-end on `test/examples/fept_tetragonal_2x2x2/input.toml`
**Script**: `julia --project test/benchmark_optimize.jl --input test/examples/fept_tetragonal_2x2x2/input.toml --with-fit --samples 5`
**Spinconfigs**: 30, salc key groups: 31, num atoms: 16.

| Section | Median | Mean ± σ | Memory | Allocs |
|---|---|---|---|---|
| `calc_∇ₑu`                            | 4.625 μs   | 8.167 μs ± 6.315 μs   | 6.97 KiB  | 236        |
| `build_design_matrix_torque` (1 cfg)  | 6.150 ms   | 6.161 ms ± 31.309 μs  | 8.37 MiB  | 318,911    |
| **`fit_sce_model` (full dataset)**    | **266.811 ms** | 266.556 ms ± 2.099 ms | 348.87 MiB | 14,331,864 |

Tag this baseline as `pre-estimator-dispatch`. Re-run after M6 with the same script and compare `fit_sce_model` median (the two helper sections are orthogonal to the refactor and just sanity checks).

## estimator-dispatch refactor — post-M6 (2026-05-13)

**Branch**: `refactor/estimator-dispatch` HEAD = `9ef0aa7` (M0–M5 applied).

**Script**: same as baseline, but **`--samples 20`** to reduce variance.

The samples=5 baseline above showed σ ≈ 0.8% of the mean inside one
run, but inter-run variation on this laptop is noticeably larger, so a
single 5-sample trial is not a reliable anchor. I re-ran the
pre-refactor baseline at samples=20 on commit `6ab646e` so the
comparison is apples-to-apples.

| `fit_sce_model` (samples=20) | Pre (`6ab646e`) | Post (`9ef0aa7`) |
|---|---|---|
| Time median            | 242.456 ms          | **240.023 ms**          |
| Time mean ± σ          | 242.250 ± 2.248 ms  | 240.288 ± 3.002 ms      |
| Memory                 | 348.87 MiB          | 348.87 MiB              |
| Allocs                 | 14,331,864          | 14,331,861              |

Within ±1% — no regression. The dispatch refactor pays no measurable
cost on the integration-style fit workload.

> Caveat: an earlier 5-sample post-refactor run on the same machine
> hit 278 ms median due to background system load. Drawing conclusions
> from 5-sample trials when σ matters is a trap — keep samples ≥ 20.

---

## B1 — CoupledBasis type parameterization (baseline)

Spec: `docs/specs/260516-coupled-basis-typeparam/`.
Branch: `refactor/coupled-basis-typeparam`.
Script: `bench/bench_b1_design_matrix.jl`.
Example: `test/examples/fept_tetragonal_2x2x2/` (num_spinconfigs=30,
num_salcs=31, num_atoms=16). 5 trials each, `@timed` + `@allocations`.

### Pre-refactor (commit `cd72f68`, branch baseline)

| function                       | min       | median    | bytes (med) | allocs (med) |
|--------------------------------|-----------|-----------|-------------|--------------|
| `build_design_matrix_energy`   | 50.9 ms   | 51.5 ms   | 101 MB      | 4,764,685    |
| `build_design_matrix_torque`   | 55.0 ms   | 55.2 ms   | 263 MB      | 9,566,937    |

### `@code_warntype` smoking gun

`design_matrix_energy_element(cbc, spin_directions, symmetry)`:

```
%46 = Base.getproperty(cbc, :coeff_tensor)::ABSTRACTARRAY
```

Downstream Any-typed locals (excerpt):

```
  Mf_size::ANY
  site_indices::ANY
  tensor_result::ANY
  mf_idx::ANY
  mf_contribution::ANY
  product_other::ANY
  m_idx_other::ANY
```

Confirms B1's premise: the bare `AbstractArray` field type forces all
downstream tensor reads to `Any`, which is responsible for the 5–10 M
allocations per call. Post-refactor target: drop both wall-time and
allocations significantly; allocation count is the cleaner signal.

### After B1 alone (Step 2)

`coeff_tensor::Array{Float64, R}` resolved; `Mf_size::Int64`,
`mf_idx::Int64` recovered. But `coeff_tensor[idx_buf..., mf_idx]` still
returns `Any` because the splat over `Vector{Int}` is not statically
resolvable. Net measured impact:

| function                       | min      | allocs (med) |
|--------------------------------|----------|--------------|
| `build_design_matrix_energy`   | 53.4 ms  | 4,701,208    |
| `build_design_matrix_torque`   | 56.1 ms  | 9,127,257    |

Wall time in noise (~+3%), allocs −1.3% / −4.6%. Step 2 was a
foundation for Step 3, not a standalone win.

### After B1 + B3 (Step 3)

Step 3 made all scratch buffers and dims compile-time-sized:

- `dims_t = ntuple(i -> 2*cbc.ls[i]+1, Val(R-1))::NTuple{R-1,Int}`
- `idx_buf, translated_atoms, atoms_sorted_buf::MVector{R-1, Int}`
- `other_dims_buf, other_sites_buf::MVector{R-2, Int}`
- `CartesianIndices` now have statically known rank

Result on the same fept fixture:

| function                       | min       | allocs (med) | speedup |
|--------------------------------|-----------|--------------|---------|
| `build_design_matrix_energy`   |   2.1 ms  |    261,625   |  ×24.2  |
| `build_design_matrix_torque`   |   6.4 ms  |  1,172,217   |  ×8.6   |

Memory: energy 101 MB → 10.1 MB, torque 263 MB → 52.6 MB. Side
effect: `make test-integration` for FeGe B20 2x2x2 dropped 26.9s →
4.2s (~6.4×) end-to-end.

Post-refactor `code_warntype` on `design_matrix_energy_element`:
`ANY count: 0` in locals; `tensor_result`, `mf_contribution`,
`product_other` all `Float64`; `idx_buf::MVector{R-1,Int64}`;
`other_site_indices::CartesianIndices{R-2, NTuple{R-2,Int64}}`.
Remaining `Union{Nothing, Tuple{Int64,Int64}}` entries are just
`iterate(::UnitRange)` return types — normal, type-stable.

---

## C — hot-path scratch workspace (`design_matrix_*_element`)

Spec: `docs/specs/260516-optimize-workspace/`.
Branch: `refactor/optimize-workspace`.

Goal: eliminate the ~262K (energy) and ~1.17M (torque) per-call heap
allocations remaining after B1+B3. Per-call profile of
`design_matrix_energy_element` showed 74-110 allocs/call (depending on
N) — dominated by per-call `Set{UInt}()` construction, per-call
`Vector{Vector{Float64}}` for spherical harmonics, and (the surprise
finding) per-(l, m) Legendre-cache allocations inside the unbuffered
`Zₗₘ_unsafe` overload.

Added `EnergyWorkspace` and `GradWorkspace` mutable structs that pool
the `Set`, the SH-value `Vector{Vector{Float64}}`, the
`atom_grad_values::Vector{SVector{3, Float64}}` (grad only), and a
`legendre_buf::Vector{Float64}` for the buffered SH overloads. Each
thread in `build_design_matrix_*` allocates one workspace inside its
`@threads` iteration.

Results on fept_tetragonal_2x2x2 (30 spinconfigs × 31 salcs):

| function                      | min       | allocs (med) | bytes (med) |
|-------------------------------|-----------|--------------|-------------|
| `build_design_matrix_energy`  | **1.5 ms**|  **9,784**   |  478 KB     |
| `build_design_matrix_torque`  | **4.1 ms**| **112,827**  |  6.95 MB    |

Versus the B1+B3 baseline (energy 2.1 ms / 262K allocs / 10.1 MB;
torque 6.4 ms / 1.17M / 52.6 MB), this is ×1.4-1.6 wall-time speedup
and ×10-27 allocation reduction. Per-call allocations dropped from
74-110 → 2 (energy) and 10-50 → 6 (grad).

Total stack vs. baseline before B1:

| function                      | pre-B1 baseline | after C |
|-------------------------------|-----------------|---------|
| `build_design_matrix_energy`  | 51.5 ms / 4.76M | 1.5 ms / 9.8K (×34) |
| `build_design_matrix_torque`  | 55.2 ms / 9.57M | 4.1 ms / 113K (×13) |

The biggest single contributor to C was the buffered SH switch —
switching `Zₗₘ_unsafe(l, m, uvec)` to `Zₗₘ_unsafe(l, m, uvec, buf)`
(same for `∂ᵢZlm_unsafe`) eliminated dozens of per-call allocations
that the workspace alone had not addressed. The buffered overload was
already in `MySphericalHarmonics.jl` (commit `8a4a17d`); only the
caller side had to wire up the workspace-owned `legendre_buf`.

## B2 — `sh_values` flatten (workspace `Vector{Vector{Float64}}` → `Vector{Float64}`)

Replaces the per-site nested layout with a single contiguous
`Vector{Float64}` plus a cumulative-offset `Vector{Int}` in
`EnergyWorkspace` / `GradWorkspace`. The kernel reads/writes
`sh_values[sh_offsets[i] + m_idx]` instead of `sh_values[i][m_idx]`.
Motivation: the original `Vector{Vector{Float64}}` had outer pointer
indirection (outer length is `R - 1`, statically known), heap-individual
inner vectors with no cache locality, and an N-iteration resize check
on every element call.

Results on fept_tetragonal_2x2x2 (30 spinconfigs × 31 salcs, 5 trials each):

| function                      | before (med allocs / time min) | after            |
|-------------------------------|--------------------------------|------------------|
| `build_design_matrix_energy`  | 9,769 / 2.3 ms                 | 9,711 / 2.3 ms   |
| `build_design_matrix_torque`  | 112,812 / 11.2 ms              | 112,782 / 11.0 ms|

Time delta sits inside BenchmarkTools noise (~±5%); allocations drop
slightly (-58 energy, -30 torque) reflecting the removed inner-vector
`push!` / `resize!` cycle in `_ensure_sh_buffer!`. The substantive
change is structural (one heap object instead of `R - 1 + 1`,
contiguous memory for the cache hierarchy) — the runtime is already
dominated by the contraction loop and SH evaluation.

Note: the comparison is against the post-workspace-pooling baseline
(`260516-optimize-workspace`), which had already removed the
per-translation allocation pattern that B2 was originally framed
around. This commit handles the residual structural concern.

## Energy-centered design matrix (spec 260518) -- `build_design_matrix_energy` shape change

Removes the bias column from `build_design_matrix_energy` /
`SCEDataset.X_E`: shape goes from `(num_spinconfigs, num_salcs + 1)` to
`(num_spinconfigs, num_salcs)` and the `design_matrix[:, 1] .= 1.0`
initialization is dropped. The bias term `j0` is recovered analytically
in `extract_j0_jphi` from the unscaled energy data, so the inner SH
contraction loop is unchanged (still iterates `i = 1:num_salcs`); only
the allocation size shrinks by one column.

Logged because `Fitting.jl` is on the hot-path list. Numerical result
is mathematically equivalent for OLS / Ridge up to FP rounding noise
(`scecoeffs.xml` diffs at the 12-15th significant digit).

Measurement (`bench/bench_b1_design_matrix.jl`, fept_tetragonal_2x2x2:
30 spinconfigs x 31 salcs x 16 atoms, 5 trials per call, 2 independent
samples reported as `a / b`):

| function                      | before (6fa751a) min          | after (5f79f84) min           | allocs (med, both)  | bytes (med, both) |
|-------------------------------|-------------------------------|-------------------------------|---------------------|-------------------|
| `build_design_matrix_energy`  | 2.0 ms / 2.1 ms               | 2.3 ms / 2.0 ms               | 9,711               | 475 KB            |
| `build_design_matrix_torque`  | 10.8 ms / 10.9 ms             | 11.1 ms / 11.1 ms             | 112,782             | 6.95 MB           |

Wall-time deltas (a few hundred us on energy, ~300 us on torque) sit
inside the run-to-run variance of single-trial `@timed` measurements;
across the two samples neither side is consistently faster. The
allocation count and total bytes are identical to the digit, as
expected -- the `(n_E, num_salcs + 1) -> (n_E, num_salcs)` matrix
shrinks by only `8 * n_E` bytes (~240 B), swamped by the SH-evaluation
allocations in the inner contraction loop. `build_design_matrix_torque`
is unchanged by this commit and serves as a noise reference.

Conclusion: performance-neutral. The refactor is shape-only on the
hot path; the substantive change lives in `assemble_weighted_problem`
(call-once-per-fit, not hot).

## Cluster-generation perf (spec 260523) -- `Cluster` construction

Two changes in `src/Clusters.jl`:

1. `irreducible_clusters` switches from an O(N_clusters^2) linear scan
   against accepted representatives (each scan-step calling
   `is_translationally_equiv_cluster`, itself iterating over all
   pure-translation symmetry ops) to a `Set{Vector{Int}}` keyed by the
   lex-minimum translation image of each cluster (`_translation_canonical_form`).
   Same equivalence relation, O(N_clusters * N_translations *
   cluster_size) instead of O(N_clusters^2 * N_translations * cluster_size).
   First-seen-wins representative is preserved, so
   `irreducible_cluster_dict` keys / order are bit-identical.
2. `set_mindist_pairs` is computed once in the `Cluster` constructor and
   threaded into `generate_clusters` as the 5th positional argument,
   eliminating a duplicated `O(N_atoms^2 * 27)` recompute.

`make test-all` (23075 tests) passes including 14367 save/load
round-trip assertions and all integration tests; numerical output is
bit-identical.

### Per-stage `Cluster` construction, `make bench-cluster` (`@timed`, sec)

Before: `bb1af45` (`bench/cluster-generation-benchmark` tip, pre-perf).
After: `da83dda` (M3 commit, all perf changes landed).

`fege_2x2x2_3body_light` (body-2 + body-3 all at 3.0 Angstrom):

| stage                  | before [s] | after [s] | speedup |
|------------------------|-----------:|----------:|--------:|
| `generate_clusters`*   |     0.0195 |    0.0028 |    7.0x |
| `irreducible_clusters` |     0.0076 |    0.0004 |   19.0x |
| `cluster_orbits`       |     0.0070 |    0.0070 |    1.0x |
| `set_mindist_pairs`    |     0.0219 |    0.0153 |    1.4x |
| **TOTAL**              | **0.0560** |**0.0237** |**2.4x** |

`fege_2x2x2_3body_fefe_open` (body-3 Fe-Fe opened, others at 3.0):

| stage                  | before [s] | after [s] | speedup  |
|------------------------|-----------:|----------:|---------:|
| `generate_clusters`*   |     0.0124 |    0.0028 |     4.4x |
| `irreducible_clusters` |     0.4284 |    0.0004 | **~1070x** |
| `cluster_orbits`       |     0.0438 |    0.0513 |     0.9x |
| `set_mindist_pairs`    |     0.0180 |    0.0106 |     1.7x |
| **TOTAL**              | **0.5027** |**0.0651** | **7.7x** |

`fege_2x2x2_3body_all_open` (all body-3 cutoffs at -1, the original
4000 s pain point reported by the user; before-figure is the user
report, not re-measured this session):

| stage                  | before [s]   | after [s] | speedup |
|------------------------|-------------:|----------:|--------:|
| `generate_clusters`*   |   (n/a)      |    0.0109 |       — |
| `irreducible_clusters` |   (dominant) |    0.0018 |       — |
| `cluster_orbits`       |   (n/a)      |    0.3042 |       — |
| `set_mindist_pairs`    |   (n/a)      |    0.0251 |       — |
| **TOTAL**              | **~4000**    |  **0.342**|**~12000x** |

*`generate_clusters` "before" includes the now-removed internal
`set_mindist_pairs` call; the post-M3 timing measures only the cluster
generation work itself.

### Bottleneck shift

`irreducible_clusters` was 85% of `fefe_open` and presumably the same
share of `all_open`. After this spec lands, `cluster_orbits` becomes
the dominant stage (89% of `all_open`) because its BFS still allocates
`sort`-fresh `Vector{Int}` per visited translation and uses a
`Set{Vector{Int}}` membership test. Deferred as a follow-up; the
remaining absolute time on `all_open` (~0.3 s) is well below the
user's pain threshold so the follow-up is not gating.

## 2026-05-23 (spec 260523-design-matrix-3body-perf)

### `build_design_matrix_energy` / `build_design_matrix_torque`, FeGe B20 2x2x2

Bench: `bench/bench_b1_3body_fege.jl --fixture light --ntrials 5`
(64 atoms, 100 spinconfigs, `JULIA_NUM_THREADS=4`, median of 5).

Before: baseline at branch tip of `main` (commit `c544b86`).
After: this spec's changes (M2 `searched_pairs` `Set` → `Vector{UInt}`,
M3 `SVector{3, Float64}` column-read hoist in both the inner
`design_matrix_energy_element` / `calc_∇ₑu!` and the outer
`build_design_matrix_torque` `dir_iatom_svec` construction).

| metric                              |   before |    after | change |
|-------------------------------------|---------:|---------:|-------:|
| `build_design_matrix_energy` (s)    |     1.32 |     1.23 |   -7 % |
| `build_design_matrix_torque` (s)    |    12.08 |    10.11 |  -16 % |
| `build_design_matrix_energy` allocs |   7.02e6 |   7.02e6 |   ±0 % |
| `build_design_matrix_torque` allocs |   3.36e8 |   3.36e8 |   ±0 % |

Allocations are unchanged: the largest remaining alloc source is the
dynamic-dispatch closure boxing on the abstract element type of
`Vector{CoupledBasis_with_coefficient}` (`Profile.Allocs` attributes
it to the captured `Symmetry` field), which this spec does not
touch. The 16 % torque speedup is real wall-time, driven by:
- M2 (`Set` → `Vector{UInt}`): ~9 % from removed hash-table overhead
  (`push!` / `in` previously consumed ~24 % of torque samples).
- M3 (`SVector` column-read hoist): ~7 % combined from removing
  `@views spin_directions[:, atom]` SubArray creation in the inner
  `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` loop (both element functions) and
  in the outer `dir_iatom_svec` construction in
  `build_design_matrix_torque`.

### Investigations explored and discarded

`Profile.Allocs` showed that the dominant remaining alloc source
(~110 M of 336 M torque allocations) is the dynamic-dispatch closure
boxing produced by `for cbc in key_group::Vector{CoupledBasis_with_coefficient}`,
attributed by the profiler to the captured
`Magesty.Symmetries.Symmetry` field. Removing this requires concretely
parameterizing the SALC element type on `R` (or a runtime
function-barrier dispatch); both are structural changes beyond this
spec and are deferred to a follow-up.

Three other rewrites were tried and rejected on this same workload:
- **Scalar `Float64` gradient accumulators** in place of
  `mf_grad_contribution` / `grad_result` `MVector{3, Float64}` and
  the broadcast `.+=`. The change was numerically bit-identical
  (per-component summation order preserved) but slightly negative
  on torque time (10.11 → 10.46 s when added on top of M2 + M3).
  The `MVector{3, Float64}` broadcast is already SROA'd into
  registers by Julia's compiler; the scalar form just removed the
  SIMD-friendly 3-wide IR pattern. Reverted.
- **Rank-specialized `if R == 2 / 3 / 4` inline contractions** of
  `coeff_tensor[idx_buf..., mf_idx]`. No measurable improvement; the
  splat indexing of an `MVector{R-1, Int}` is already compiler-elided
  to a static `getindex`.
- **Generic `CartesianIndices(dims_t)` contraction with
  workspace-backed `Vector{Int}` scratch** (no per-call `MVector`).
  Allocs dropped 67 % but torque wall time regressed 8 % because the
  gradient version had to recompute `sh_prod` `2*l_atom + 1` more
  times than the original splat form. The remaining 112 M allocations
  were all dispatch-boxing, not in our hot path. Reverted; the
  simpler splat-based code in place after the revert is what produced
  the 16 % speedup above.

## `cluster_orbits` BFS allocation cleanup -- 3x3x3 FeGe pain point

The BFS inner body in `cluster_orbits` (`src/Clusters.jl`) used to
allocate two fresh `Vector{Int}` per `(symop, sym_tran)` iteration (via
list-comprehension), `sort` them, and run set-membership lookup. On a
real 3x3x3 FeGe B20 (216 atoms, body-2 and body-3 cutoffs all `-1`) the
stage allocated ~88 GiB and dominated `Cluster` construction at 98.7 %.
Refactor: two scratch `Vector{Int}` buffers of size `body` are
pre-allocated per outer body loop; the per-iteration shifts /
translations / sorts run in-place; a `copy(buf)` is taken only when a
candidate actually inserts. Algorithm and insertion order unchanged.

Bench: `bench/benchmark_cluster.jl`, single thread, `@timed`, one warmup
pass, one measured pass. Same machine session for before / after.

`fege_2x2x2_3body_all_open` (64-atom, body-3 cutoffs all `-1`):

| stage                  | before [s] | after [s] | alloc before [MiB] | alloc after [MiB] |
|------------------------|-----------:|----------:|-------------------:|------------------:|
| `set_mindist_pairs`    |     0.0251 |    0.0098 |               --   |              67.4 |
| `generate_clusters`    |     0.0109 |    0.0106 |               --   |              57.4 |
| `irreducible_clusters` |     0.0018 |    0.0018 |               --   |               4.1 |
| `cluster_orbits`       |     0.3042 |    0.0744 |               --   |               1.6 |
| **TOTAL**              | **0.342**  |**0.097**  |                    |                   |

`fege_3x3x3` real workload (216-atom B20, body-2 and body-3 cutoffs
all `-1`, the input that triggered the investigation):

| stage                  | before [s] | after [s] | alloc before [MiB] | alloc after [MiB] | speedup |
|------------------------|-----------:|----------:|-------------------:|------------------:|--------:|
| `set_mindist_pairs`    |     0.233  |    0.265  |              766.8 |             766.8 |    0.9x |
| `generate_clusters`    |     0.103  |    0.156  |              511.5 |             511.5 |    0.7x |
| `irreducible_clusters` |     0.059  |    0.054  |               41.8 |              41.8 |    1.1x |
| `cluster_orbits`       |    29.677  |   12.816  |          88 399.0  |              14.3 | **2.3x** |
| **TOTAL**              |**30.07**   |**13.29**  |                    |                   | **2.3x** |

Body-3 cluster counts on 3x3x3: raw 95832, irreducible 31944, orbits
2678. `cluster_orbits` allocation count dropped 6190x.
`set_mindist_pairs` / `generate_clusters` wall-times are unchanged
within run-to-run noise (the refactor does not touch them; the
small jitter above is single-pass `@timed` variance).

`cluster_orbits` is still 96 % of the post-refactor 3x3x3 total
(12.8 s absolute). The remaining cost is the BFS itself --
inherently sequential per orbit (frontier dependency), with shared
`processed_clusters` state across orbits. Parallel BFS (Tier 2)
would need a thread-safe claim map and lex-min orbit re-numbering to
preserve deterministic `orbit_index`. Deferred pending a dedicated
spec.

`make test-unit` (22136 passed), `make test-integration` (883
passed), `make test-jet`, `make test-aqua` all pass; numerical
output is bit-identical (same insertion order; only allocation
shape changed).

## `cluster_orbits` parallel BFS via union-find -- 3x3x3 FeGe

Building on the previous buffer-reuse cleanup. The per-orbit BFS was
inherently sequential (frontier dependency, shared `processed_clusters`
across orbits). Refactor: rewrite `cluster_orbits` (`src/Clusters.jl`)
into three phases.

- Phase 1 (parallel via `Threads.@threads`): per cluster, enumerate
  the *set* of in-set neighbor indices reachable in one (spatial
  symmetry op, translation) step. Each iteration reads `clusters_sorted`,
  `cluster_to_idx`, and `symmetry.{symdata, symnum_translation, map_sym,
  map_sym_inv}` -- all immutable in this phase -- and writes only its
  own slot of `neighbors_per_cluster`. Thread-local
  `buf_shifted` / `buf_translated` / `local_set`.
- Phase 2 (serial): union-find merge of the edges discovered in
  Phase 1. Path halving + union by size; negligible relative to Phase 1.
- Phase 3 (serial): group cluster indices by root, sort each orbit's
  members lex, sort orbits by lex-min member, assign `orbit_index` 1..N
  in that order. Deterministic and independent of thread scheduling.

**Output contract**: orbit *sets* (the partition) and the
`orbit_index -> orbit set` mapping are bit-identical to the previous
sequential BFS. The *order of cluster members within each orbit's
`Vector{Vector{Int}}`* changed from BFS visit order to lex-sorted.
No downstream consumer reads intra-orbit member order: `SALCBases.jl`
sorts atom lists internally before comparing; tests assert counts and
Set semantics. Verified via `make test-all` at 1 and 8 threads (23019
passed each).

Bench: `bench/benchmark_cluster.jl`, `@timed`, one warmup, one measured
pass per stage. Same machine session.

`fege_3x3x3` real workload (216-atom B20, body-2 and body-3 cutoffs
all `-1`):

| `cluster_orbits` [s] |  baseline | +Tier1 (buffer) | +Tier2 1T | +Tier2 4T | +Tier2 8T |
|----------------------|----------:|----------------:|----------:|----------:|----------:|
| time                 |     29.68 |          12.82  |     13.24 |      3.87 |      2.53 |
| total `Cluster` [s]  |    30.07  |          13.29  |     13.71 |      4.28 |      2.92 |
| alloc [MiB]          | 88 399.0  |           14.3  |      40.6 |      40.6 |      40.6 |
| speedup vs baseline  |     1.0x  |           2.31x |     2.24x |     7.67x | **11.73x** |

Parallel efficiency at 8 threads on the cluster_orbits stage:
13.24 / 2.53 / 8 = 65 %. Phase 1 (`Dict.get` lookups on
`cluster_to_idx`) is memory-bound and dominates the wall time;
Phase 2 (union-find) and Phase 3 (sort + emit) are serial but
small. Hoisting the per-iteration `Vector{Int}` / `Set{Int}`
scratch to `threadid()`-indexed thread-local slots with `:static`
scheduling was measured and rejected: it cut allocations from
40 MiB to 13 MiB but regressed 8-thread wall time by 10-15 %
(`empty!` plus the `bufs[tid]` indirection cost more than the
saved allocator work; the GC is not the bottleneck at this
allocation rate). The `Set{Int}` size hint
(`sizehint!(local_set, symmetry.ntran)`) is the only retained
allocation-side change -- it lets the Set skip a mid-loop rehash
without restructuring the loop.

`fege_2x2x2_3body_all_open` regression check (64-atom; small absolute
times; orbits 5 + 268):

| stage                  | +Tier2 8T [s] | alloc [MiB] |
|------------------------|--------------:|------------:|
| `set_mindist_pairs`    |        0.0107 |         67.4 |
| `generate_clusters`    |        0.0196 |         57.4 |
| `irreducible_clusters` |        0.0020 |          4.1 |
| `cluster_orbits`       |        0.0172 |          4.7 |
| **TOTAL**              |    **0.0495** |              |

On the small case `generate_clusters` (sequential) is now the
dominant stage; `cluster_orbits` has dropped to 35 % and is no longer
gating. No regression on the regression fixture.

`Cluster.cluster_orbits_dict` allocation under the new code is
~48 MiB on 3x3x3 (independent of thread count): roughly the
`Set{Int}` neighbor buffers + the `Dict{Vector{Int}, Int}` index map
+ the final orbit `Vector{Vector{Int}}` shells. This is 1840x lower
than the original 88 GiB, with most of the residual coming from the
`Set{Int}` per-thread neighbor accumulation in Phase 1.
