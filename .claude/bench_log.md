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
