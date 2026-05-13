# Optimize.jl の estimator dispatch リファクタリング

**Status**: **完了** (2026-05-13、branch `refactor/estimator-dispatch`)

Spec: `docs/specs/260513-estimator-dispatch/`。

## ハイライト

- 3 層分離（`assemble_weighted_problem` → `solve_coefficients` → `extract_j0_jphi`）と多重ディスパッチに移行。
- `ElasticNet`（実体は pure ridge）→ `Ridge(lambda)` に **rename**（breaking）。`ElasticNet` 名前を将来の真の mixed-norm estimator のために空けた。
- `[regression].alpha` は TOML 互換のため field/parse は残し、`alpha != 0` でワンショットの deprecation `@warn` を発行（`maxlog = 1`）。
- Bench: `fit_sce_model` (`fept_tetragonal_2x2x2`, samples=20) で +/-1%（`.claude/bench_log.md` 参照）。

以下は当初の設計メモ（参考用）。

---

**対象**: `src/Optimize.jl` の回帰手法ディスパッチと `elastic_net_regression`。

**動機**: 将来的に Lasso / Ridge / Bayesian / NNLS / 制約付き回帰などを追加する可能性がある。現状は型階層（`AbstractEstimator` / `OLS` / `ElasticNet`）は存在するが、`_fit_sce_model_internal` で `isa` 分岐しており、多重ディスパッチの利点を活かしていない。さらに `elastic_net_regression` が「重み付き問題組み立て」「求解」「j0 抽出」を 1 関数に混在させているため、新しい回帰を足すと組み立て・後処理がコピペで増殖する。

## 現状の構造的弱点

`elastic_net_regression`（L938–998）が独立した 3 責務を抱えている:

1. **問題組み立て**（estimator 非依存）— energy/torque を √weight でスケール、縦連結、bias 列整合。
2. **求解**（estimator 依存）— `X \ y` または `ridge(X, y, λ; bias=false)`。
3. **後処理**（estimator 非依存）— `j0 = mean(y_E - X_E[:,2:end] * jphi)` で bias 抽出。

`_fit_sce_model_internal` の `isa` 分岐（L837–857）は型を持ちながら dispatch しない「両方の悪いとこ取り」状態。

## 提案: 3 層分離 + 多重ディスパッチ

```julia
abstract type AbstractEstimator end

struct OLS        <: AbstractEstimator end
struct Ridge      <: AbstractEstimator; lambda::Float64 end
struct Lasso      <: AbstractEstimator; lambda::Float64 end
struct ElasticNet <: AbstractEstimator; alpha::Float64; lambda::Float64 end
# struct NNLS, ConstrainedOLS, ElasticNetCV, ...

# 第1層: 問題組み立て（estimator 非依存）
assemble_weighted_problem(Xe, Xt, ye, yt, weight) -> (X, y, bias_col)

# 第2層: 求解（estimator ごとに dispatch、新しい回帰の追加点はここだけ）
solve_coefficients(::OLS,        X, y; bias_col) = X \ y
solve_coefficients(e::Ridge,     X, y; bias_col) = ridge_solve(X, y, e.lambda; bias_col)
solve_coefficients(e::Lasso,     X, y; bias_col) = lasso_solve(X, y, e.lambda; bias_col)
solve_coefficients(e::ElasticNet,X, y; bias_col) = enet_solve(X, y, e.alpha, e.lambda; bias_col)

# 第3層: 係数抽出（estimator 非依存）
extract_j0_jphi(j_values, Xe, ye) -> (j0, jphi)
```

これで `_fit_sce_model_internal` は 3 行に縮む:

```julia
X, y, bias_col = assemble_weighted_problem(Xe, Xt, ye, yt, weight)
j_values        = solve_coefficients(estimator, X, y; bias_col)
return extract_j0_jphi(j_values, Xe, ye)
```

**新しい回帰の追加 = 新しい struct + `solve_coefficients` の 1 メソッドだけ**（open/closed 原則）。`isa` 分岐の `else throw(ArgumentError(...))` も不要（未定義型は `MethodError` で自動失敗）。

## 設計上の判断ポイント

- **bias 列の扱い**: 「bias を正則化から除外」は estimator 共通要件。`bias_col` を引数で渡すことで、各 solver が自分の方法で除外できる（自前実装なら `lambda_vec[bias_col]=0`、GLMNet.jl 系なら `intercept=true`）。今のうちにこのインターフェースを引き出しておく。
- **反復解法のオプション**（`max_iter`, `tol`, `warm_start`）: 各 estimator struct の field として持たせる。Holy traits は estimator が 10 を超えるまで導入しない。
- **CV / λ-path**: `ElasticNetCV` のような**別の estimator 型**として `solve_coefficients` を定義する。`ElasticNet` の中に CV ロジックを混ぜない。
- **外部パッケージ依存**（Lasso.jl, GLMNet.jl 等）: `solve_coefficients` の内部に閉じ込め、上位 API には漏らさない。

## 付随で整理したい点

- `ElasticNet.alpha` が "currently unused"（L46, L931）。L1 を入れる予定がないなら型名を `Ridge` に変えた方が誤解がない。L1 を入れる予定なら docstring に「現状 α≠0 でも無視される」を明記。
- `Optimizer` 構造体コンストラクタ（L94–104）で `alpha`, `lambda` を位置引数で取りつつ `estimator::AbstractEstimator = ElasticNet(alpha=alpha, lambda=lambda)` をデフォルトにしている設計は曖昧。estimator を渡したら位置引数が無視されるのか上書きされるのか呼び出し側から読み取れない。リファクタと合わせて `Optimizer(..., estimator, weight)` に寄せる。

## やらないこと（現時点）

- StatsAPI.jl / MLJ.jl 互換化（外部公開時に検討、今は内部 API で十分）。
- Holy traits（`CapabilityTrait` 系）— 4–5 種類なら多重 dispatch で捌ける。
- `AbstractEstimator{T}` のような抽象パラメータ化（YAGNI）。

## 進め方

API シグネチャに影響する中規模リファクタなので、着手前に `docs/specs/[YYMMDD]-estimator-dispatch/` を切って requirements.md / design.md / tasklist.md で合意する。実装単独で着手しない。

## 連動箇所

- `_fit_sce_model_internal`, `fit_sce_model_ols`, `fit_sce_model_elastic_net`, `elastic_net_regression`（L829–998）
- `Optimizer` 構造体の外側 constructor（L94–283）
- `fit_sce_model`（L678–790）の `estimator` 引数まわり
- `export` 一覧（L24）
