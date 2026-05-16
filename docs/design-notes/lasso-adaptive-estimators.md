# LASSO / Adaptive LASSO / Adaptive Ridge estimator の導入

**Status**: 未着手 (2026-05-16)

GLMNet.jl を採用し、`Optimize.jl` の estimator dispatch ファミリに L1 系・Adaptive 系を追加する設計メモ。spec 化前の合意ドラフト。

## 背景

現状 `src/Optimize.jl` の SCE 係数推定は `OLS`（`X \ y`）と `Ridge`（`MultivariateStats.ridge` + per-column penalty で bias 列を非ペナルティ化）の 2 つのみ。スパースな有効スピンモデルを得るために L1 正則化系（LASSO, Adaptive LASSO）が欲しく、また L0 近似として Adaptive Ridge（Frommlet–Nuel 2016 の反復重み付け方式）も導入したい。

## ライブラリ選定

`MultivariateStats.jl` には LASSO は無い。Julia エコシステムの候補比較:

| 観点 | Lasso.jl | GLMNet.jl |
|---|---|---|
| 実装 | 純 Julia (coordinate descent) | Friedman/Hastie/Tibshirani の Fortran `glmnet` ラッパ（reference 実装） |
| 信頼性 | JuliaStats 公式・活発にメンテ | R glmnet・scikit-learn とビット同等。30 年の実績 |
| Julia 依存 | GLM, StatsModels, StatsBase, Distributions 等やや多い | Distributions, SparseArrays + GLMNet_jll のみ |
| JLL 依存 | 無し | あり（本 project は Spglib_jll 等で既存パターン、追加負担なし） |
| `penalty_factor`（bias 除外・Adaptive 重み） | 対応 | 対応 |
| `alpha` で LASSO / Ridge / Elastic Net 統一 | 一部 | 完全対応 |

**採用: GLMNet.jl**。理由は (a) reference 実装としての数値再現性（CLAUDE.md の「数値的な正確さ・再現性を優先」と合致）、(b) Julia 依存が最小、(c) `alpha` と `penalty_factor` だけで LASSO / Adaptive LASSO / Adaptive Ridge / Elastic Net を一つの実装で網羅できる構造的単純さ。

既存 `Ridge` は MultivariateStats のまま変更しない（リグレッションリスク回避）。新規 estimator のみ GLMNet 経由。

## 規約（既存 Ridge と揃える）

| 項目 | 方針 |
|---|---|
| bias 列 | `penalty_factor[bias_col] = 0` で非ペナルティ化 |
| `j0` 抽出 | 既存の `extract_j0_jphi` をそのまま再利用 |
| standardize | デフォルト OFF（SCE スケール因子で列が既に物理的に意味づけられている）。`standardize::Bool` kwarg で ON 可 |
| λ 選択 | 単一 λ をユーザ指定（Ridge と一貫）。CV 自動選択は別 spec で扱う |
| augmented system | `assemble_weighted_problem` の出力 `(X, y, bias_col)` をそのまま渡す |

## estimator マッピング

| Estimator | `alpha` | `penalty_factor[j ≠ bias]` | 補足 |
|---|---|---|---|
| Lasso | 1.0 | 1 | — |
| AdaptiveLasso | 1.0 | `1/|β̂_init,j|^γ` | β̂_init は OLS or Ridge（Zou 2006） |
| AdaptiveRidge (`mode = :oneshot`) | 0.0 | `1/|β̂_init,j|^γ` | 一段階重み付け |
| AdaptiveRidge (`mode = :iterative`) | 0.0 | `1/(β_j² + ε)` を反復更新 | Frommlet–Nuel 2016。各反復で GLMNet を再 fit して L0 近似稀疏化 |

`:oneshot` / `:iterative` は動作記述命名（人名 symbol は読み手が論文を知る前提になるため避けた）。

## API 草案

```julia
struct Lasso <: AbstractEstimator
    lambda::Float64
    standardize::Bool
end

struct AdaptiveLasso <: AbstractEstimator
    lambda::Float64
    gamma::Float64           # 典型 1.0
    init::Symbol             # :ols | :ridge
    init_lambda::Float64     # init = :ridge のみ使用
    standardize::Bool
end

struct AdaptiveRidge <: AbstractEstimator
    lambda::Float64
    mode::Symbol             # :oneshot | :iterative
    gamma::Float64           # :oneshot のみ
    init::Symbol             # :oneshot のみ (:ols | :ridge)
    init_lambda::Float64
    epsilon::Float64         # :iterative のみ（β_j² への加算）
    max_iter::Int            # :iterative のみ
    tol::Float64             # :iterative のみ（係数 L∞ 変化）
    standardize::Bool
end
```

`solve_coefficients` 各メソッドは共通ヘルパー `_glmnet_solve(X, y, alpha, penalty_factor; ...)` を呼ぶ。`_initial_estimate` / `_adaptive_penalty_factor` / `_iterative_penalty_factor` を内部ヘルパーとして用意。

`Lasso` / `AdaptiveLasso` / `AdaptiveRidge` を export し `Magesty.jl` で re-export。

## テスト計画

1. **λ → 0 で OLS と一致**（要件）: 小規模 well-conditioned `X, y` で `Lasso(λ=1e-10)` / `AdaptiveLasso(λ=1e-10)` / `AdaptiveRidge(λ=1e-10, mode=:oneshot)` / `AdaptiveRidge(λ=1e-10, mode=:iterative)` の係数が `OLS()` と `atol` 内で一致。bias 係数（j0）も含めて比較。
2. **bias 列が penalty から除外される**: penalty を極端に大きくしても bias_col 係数は `mean(y)` 相当（他がゼロに落ちる極限）。
3. **Adaptive LASSO の oracle 性**: スパース真解の合成データでゼロ係数が実際にゼロに落ちる。
4. **AdaptiveRidge `:iterative` の収束**: `max_iter` 未満で `tol` を満たす。
5. **GLMNet と MultivariateStats.ridge の整合**: `alpha=0, penalty_factor=unit` の GLMNet 解が既存 Ridge 解と一致（数値再現性）。
6. `make test-all` `make test-jet` `make test-aqua` を通す。

## 修正対象

- `Project.toml` — `[deps]` `[compat]` に `GLMNet`
- `src/Optimize.jl` — `using GLMNet`、3 つの型・メソッド・ヘルパー・export
- `src/Magesty.jl` — re-export
- `test/component_test/test_optimize_estimators.jl`（新規 or 既存に追記）

## 進め方

中規模機能追加 + 設計判断複数 + 新外部依存 → spec 化が妥当。合意後に `docs/specs/260516-lasso-adaptive-estimators/`（requirements / design / tasklist の 3 ファイル）を切り、ブランチ `refactor/lasso-adaptive-estimators` で実装する（[[feedback_branch_for_medium_refactor]]）。

## オープン項目

- `standardize` デフォルト OFF でよいか（Ridge と一貫させる方針）。
- `AdaptiveLasso` / `AdaptiveRidge(:oneshot)` の `init` デフォルト `:ols`、`γ` デフォルト `1.0`。
- `AdaptiveRidge(:iterative)` のデフォルト `epsilon=1e-8`, `max_iter=50`, `tol=1e-6`。
- 既存 `Ridge` を将来 GLMNet に寄せて API を一本化するかは別 spec で判断。
