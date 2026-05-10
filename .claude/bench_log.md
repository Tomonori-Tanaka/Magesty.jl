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
