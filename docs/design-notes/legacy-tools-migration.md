# Legacy tools migration

Status: 保留中（作成 2026-05-16）。

`tools/` 配下の以下のスクリプトは Step 7 (SCE 公開 API 破壊コミット) 以前の
シンボル (`System`, `SpinCluster`, `fit_sce_model`,
`_fit_sce_model_internal`, `get_j0_jphi`, `calc_energy` 等) を参照しており、
現状の `Magesty.jl` ロード時点で `UndefVarError` になる。

## 対象ファイル

| ファイル | 行数 | 用途 |
|---|---|---|
| `tools/check_convergence_embset.jl` | 394 | EMBSET の train 数を `--n-range` で振って jphi / energy 予測誤差の収束を TSV/PNG 出力する CLI |
| `tools/convert2tensor.jl` | 754 | SCE 係数を交換テンソル (`{J_ij^αβ}`) 形式に変換するモジュール |
| `tools/micromagnetics.jl` | 111 | SCE からマイクロマグネティクスモデルパラメータを導出する CLI |
| `tools/plot_jphi_cluster_distance.jl` | 229 | `jphi.xml` を読んで SALC 係数 vs クラスタ距離を可視化する CLI |
| `tools/personal/*` | (gitignored) | 個別ユーザー用 — 移行はユーザー判断 |

## 移行パターン

旧 → 新の対応関係（Step 7 で確立）:

| 旧 | 新 |
|---|---|
| `System(input; verbosity = false)` | `SCEBasis(input; verbosity = false)` |
| `SpinCluster(system, input)` または `SpinCluster(system, embset_path)` | `SCEDataset(basis, ...)` |
| `fit_sce_model(spincluster, ...)` / `_fit_sce_model_internal` | `fit(SCEFit, dataset, OLS(); torque_weight = ...)` |
| `get_j0_jphi(fitted_obj)` | `intercept(fit)`, `coef(fit)` |
| `calc_energy(model, spin_config)` | `predict_energy(model, spin_config)` |
| `system.jld2` でキャッシュしたパターン | `SCEBasis(toml)` を毎回構築、または `save(basis, "basis.xml")` / `load(SCEBasis, "basis.xml")` |

## 方針

- 「使いたくなった人が」その時点で移行する。スケジュールしない。
- `tools/test/` には旧 API テストが残っている可能性があるので、各 CLI を移行する際は対応するテストも更新する。
- 移行と同時に `docs/src/tools.md` の該当節を更新する（旧 `System` / `system.jld2` 言及が残存中）。

## 削除済みの関連スクリプト（参考）

2026-05-16 にまとめて削除済み（`b163fb3 chore(bench): delete legacy benchmark and profile scripts` 参照）:

- `tools/CrossValidation.jl` — LOOCV モジュール、外部参照なし。
- `test/benchmark_optimize.jl`, `test/benchmark_optimize_hotspots.jl`, `test/benchmark_optimize_sphericart.jl`, `test/benchmark_threads.jl`, `test/profile_run.jl` — Makefile target / TEST_MODE 分岐ごと削除。

## 参考

- `docs/design-notes/post-step7-cleanup.md` — Step 7 直後の追走項目全体。Out-of-scope 節で同じ tools/ に言及している。
- `docs/specs/260514-sce-public-api/` — Step 7 の元 spec（旧 → 新 API の確定経緯）。
