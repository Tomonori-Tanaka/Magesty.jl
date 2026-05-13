# SCE 公開 API: 4 型構成 + StatsAPI 化

**Status**: 未着手（2026-05-13 提案、2026-05-14 拡張版に更新）

**結論案**: 4 つの型 `SCEProblem` / `SCEDataset` / `SCEFit` / `SCEModel` に整理し、StatsAPI.jl の語彙（`fit` / `coef` / `predict` 等）に合わせる。破壊的変更を伴うため pre-release（`0.1.0-DEV`）のうちに整理する。後方互換性より理想構成を優先。

## 背景

### 現状の公開型と問題点

| 段階 | 現名 | 中身 | 問題 |
|---|---|---|---|
| ① セットアップ | `System` | `Structure`+`Symmetry`+`Cluster`+`BasisSet` | `System` は generic すぎて `ModelingToolkit.System` 等と衝突しうる |
| ② フィット | `SpinCluster` | ① + `Optimizer`（design matrix + 訓練配置 + 係数 + メトリクス） | 内部に `Cluster` 型があるため "spin の cluster" と読まれて混乱。design matrix と係数を一体に抱えていて再利用が効きにくい |
| ③ 予測モデル | `SCEModel` | `j0`, `jphi`, `basisset`, `symmetry`, `num_atoms` | 妥当な名前。据置 |

### 想定する使い方

加えて以下の使い方が想定される:

- **weight を変えて何度もフィット**（energy vs torque tradeoff の sweep）。weight 込みの X を保持すると毎回再計算が必要になり無駄。
- **訓練データのサブセットでフィット**（train/test split, k-fold CV, ブートストラップ）。design matrix を作り直さずスライスで済ませたい。
- **複数の EMBSET を結合**（強磁性ベースと反強磁性ベースを別々に作って後で `vcat`）。
- **高スループット計算**（組成スイープ・格子定数スキャン）。`for compound in compounds` のスクリプティングが自然に書けることを設計の検証基準とする。

これらを自然に表現するには「design matrix を抱える中間型」が要る。

### 設計原則

1. **Julia API を一級市民**: 物質定義もアルゴリズムも Julia コンストラクタで完結できる。TOML はそれを呼ぶ薄いラッパー（コンビニエンス層）。
2. **StatsAPI / GLM.jl の語彙に合わせる**: Julia の ML/統計エコシステムの常識に乗せる。
3. **責務を 4 層に分離**: 物質+基底（`SCEProblem`）／データ＋design matrix（`SCEDataset`）／フィット結果（`SCEFit`）／予測器（`SCEModel`）を独立に構築・再利用できる。
4. **高スループットを一級でサポート**: parameter sweep が自然に書けることを設計の検証基準とする。

### Julia 生態系の規範

主要な回帰/最適化パッケージ:

- GLM.jl: `lm(formula, data)` → `LinearModel`、`fit(LinearModel, X, y)`
- MixedModels.jl: `fit(MixedModel, formula, data)` → 結果オブジェクト
- LsqFit.jl: `curve_fit(model, xdata, ydata, p0)` → `LsqFitResult`
- Optim.jl: `optimize(f, x0, method)` → `OptimizationResults`
- DifferentialEquations.jl: `solve(problem, alg)` → `ODESolution`

StatsAPI.jl で `fit` / `coef` / `predict` / `residuals` 等が抽象動詞として定義されており、軽量依存（純粋関数の抽象定義のみ）で乗りやすい。

## 4 型の役割分担

| 型 | 役割 | 中身（概念） | 重さ |
|---|---|---|---|
| `SCEProblem` | 物質+基底の静的セットアップ | `Structure`+`Symmetry`+`Cluster`+`BasisSet` | 軽量 |
| `SCEDataset` | 訓練データを X 化した状態 | `problem`+`X_E`+`X_T`+`y_E`+`y_T`+`spinconfigs`（unweighted）| 重い・再利用可能 |
| `SCEFit` | フィット結果（診断用フル装備） | `dataset` 参照+係数+メトリクス+残差+使った estimator/weight | 中 |
| `SCEModel` | 軽量予測器（シリアライズ用） | `j0`+`jphi`+`basisset`+`symmetry`+`num_atoms` | 軽量 |

**重要な設計判断**: `SCEDataset` は **weight を掛ける前の** `X_E` / `X_T` を保持する。weight は `fit` 時のパラメタとして渡し、weight sweep で X を作り直さずに済むようにする。

## 提案する使用例スクリプト

### 例 1: 基本フロー

```julia
using Magesty

problem = SCEProblem("input.toml")                              # 物質+基底（軽量）
spinconfigs = read_embset("EMBSET.dat")
dataset = SCEDataset(problem, spinconfigs)                      # X 構築（重い、再利用可能）

# OLS
fit_ols = fit(SCEFit, dataset, OLS())
println("RMSE energy: ", metrics(fit_ols)[:rmse_energy] * 1000, " meV")

# Ridge with weight = 0.5（同じ dataset を再利用、X 再構築なし）
fit_ridge = fit(SCEFit, dataset, Ridge(lambda=1e-4); weight=0.5)

# 予測器に変換して保存
model = SCEModel(fit_ridge)
save(model, "model.xml")
```

### 例 2: weight sweep（design matrix を 1 回構築して使い回す）

```julia
problem = SCEProblem("input.toml")
dataset = SCEDataset(problem, read_embset("EMBSET.dat"))         # X 構築は 1 回のみ

weights = 0.0:0.1:1.0
fits = [fit(SCEFit, dataset, Ridge(lambda=1e-4); weight=w) for w in weights]

rmse_e = [metrics(f)[:rmse_energy] for f in fits]
rmse_t = [metrics(f)[:rmse_torque] for f in fits]
```

### 例 3: train/test split（dataset のスライス）

```julia
dataset = SCEDataset(problem, read_embset("EMBSET.dat"))

n_train = round(Int, 0.8 * length(dataset))
train, test = dataset[1:n_train], dataset[n_train+1:end]

fit_result = fit(SCEFit, train, Ridge(lambda=1e-4); weight=0.5)
model = SCEModel(fit_result)

# テストセットで評価
test_rmse = sqrt(mean((predict(model, sc.spin_directions) - sc.energy)^2
                       for sc in test.spinconfigs))
```

### 例 4: k-fold CV（自作ヘルパーで kfold 関数を作る想定）

```julia
folds = kfold(dataset, 5)   # Vector{Tuple{SCEDataset, SCEDataset}}
cv_rmse = mean(folds) do (train, val)
    fit_result = fit(SCEFit, train, Ridge(lambda=1e-4); weight=0.5)
    model = SCEModel(fit_result)
    sqrt(mean((predict(model, sc.spin_directions) - sc.energy)^2
              for sc in val.spinconfigs))
end
```

### 例 5: 複数 EMBSET の結合（FM + AFM）

```julia
problem = SCEProblem("input.toml")

fm_data  = SCEDataset(problem, read_embset("EMBSET_fm.dat"))
afm_data = SCEDataset(problem, read_embset("EMBSET_afm.dat"))

# 同じ problem に紐づくことを runtime check（basisset/num_atoms 一致）
combined = vcat(fm_data, afm_data)

fit_result = fit(SCEFit, combined, Ridge(lambda=1e-4); weight=0.5)
```

### 例 6: estimator 比較（同じ dataset を使い回す）

```julia
dataset = SCEDataset(problem, read_embset("EMBSET.dat"))

estimators = [OLS(), Ridge(lambda=1e-4), Lasso(lambda=1e-3),
              ElasticNet(alpha=0.5, lambda=1e-4)]
fits = [fit(SCEFit, dataset, est; weight=0.5) for est in estimators]

for (est, f) in zip(estimators, fits)
    println(typeof(est), ": rmse_E = ", metrics(f)[:rmse_energy] * 1000, " meV")
end
```

### 例 7: 高スループット — 格子定数スイープ

```julia
results = map(3.0:0.05:3.5) do a
    crystal = Crystal(; lattice = a * I(3), kinds=["Fe","Fe"],
                       positions=[SVector(0,0,0), SVector(0.5,0.5,0.5)])
    problem = SCEProblem(crystal; nbody=2, lmax=Dict("Fe"=>2))
    dataset = SCEDataset(problem, read_embset("EMBSET.dat"))
    fit(SCEFit, dataset, Ridge(lambda=1e-4); weight=0.5)
end
```

### 例 8: 基底を XML からロードして再フィット

```julia
problem = SCEProblem("input.toml"; basis_from = "scecoeffs.xml")
dataset = SCEDataset(problem, read_embset("new_EMBSET.dat"))
fit_result = fit(SCEFit, dataset, Ridge(lambda=1e-3); weight=0.5)
save(SCEModel(fit_result), "new_results.xml")
```

### 例 9: nbody / lmax 収束スキャン

```julia
crystal = Crystal(; lattice=a*I(3), kinds=["Fe","Fe"],
                   positions=[SVector(0,0,0), SVector(0.5,0.5,0.5)])
data = read_embset("EMBSET.dat")

convergence_runs = [
    (nbody=2, lmax=Dict("Fe"=>1)),
    (nbody=2, lmax=Dict("Fe"=>2)),
    (nbody=2, lmax=Dict("Fe"=>3)),
    (nbody=3, lmax=Dict("Fe"=>2)),
    (nbody=3, lmax=Dict("Fe"=>3)),
]
results = map(convergence_runs) do params
    problem = SCEProblem(crystal; params..., tolerance_sym=1e-3)
    dataset = SCEDataset(problem, data)
    fit(SCEFit, dataset, Ridge(lambda=0.1))
end
rmses = [metrics(f)[:rmse_energy] for f in results]
```

これらは現状の input.toml フローでは TOML ファイルを書き換えて Julia を再起動する必要があるが、本 API では純粋な Julia ループで完結する。

## スライス API のセマンティクス

- `length(dataset)` → spinconfig 数（n_configs）
- `dataset[5]` → 5 番目の spinconfig 1 個を含む `SCEDataset`（型保存的）
- `dataset[1:80]` → 最初 80 spinconfigs を含む `SCEDataset`
- `dataset[mask::BitVector]` → mask で選択した spinconfigs を含む `SCEDataset`
- `vcat(d1, d2)` → 結合した `SCEDataset`（同じ `problem` 必須、runtime check）

`X_E` の行スライスと `X_T` の対応ブロック（`3·n_atoms` 行）スライスを内部で同期する。

## StatsAPI 準拠の動詞

```julia
coef(fit_or_model)        # → jphi::Vector{Float64}
intercept(fit_or_model)   # → j0::Float64
predict(model, spin_dir)  # → energy::Float64
residuals(fit)            # → Vector{Float64}
metrics(fit)              # → Dict{Symbol,Any}
dataset(fit)              # → SCEDataset
```

`predict_torque(model, spin_dir)` は StatsAPI に対応動詞がないので Magesty 独自。

## 拡張性: Cross-Validation / Bayesian への対応

「将来 Lasso 以外も入れる」想定で本設計を試したところ、**Cross-Validation と Bayesian 線形回帰の両方ともこのアーキテクチャに自然に乗る**。ただし dispatch レベルと result type に関する設計判断が必要なので、明記しておく。

### Cross-Validation の組み込み方

CV は **`solve_coefficients` レベルではなく `fit` レベルの override** にする。理由:

- `solve_coefficients(est, X, y; bias_col)` は assembled な X, y（energy + torque を縦連結したもの）で動く。
- CV は spinconfig 単位で fold 分割が必要だが、assembled X, y では 1 spinconfig が `1 + 3·num_atoms` 行に展開済み。行レベルでランダム分割すると意味が壊れる。
- 本 4 型構成では `SCEDataset` がスライス可能なので、CV ヘルパは dataset 上で fold を作る形になる。

そこで CV は wrapper estimator として実装し、`fit` 側で dispatch する:

```julia
struct CV{E<:AbstractEstimator} <: AbstractEstimator
    base::E
    k::Int
    seed::Int
end

function fit(::Type{SCEFit}, dataset::SCEDataset, cv::CV; weight)
    folds = kfold(dataset, cv.k; seed=cv.seed)
    scores = map(folds) do (train, val)
        f = fit(SCEFit, train, cv.base; weight)
        evaluate(SCEModel(f), val)
    end
    # 必要ならハイパラ探索を行い、全データで refit して返す
    return fit(SCEFit, dataset, cv.base; weight)
end
```

ユーザーコードは `fit(SCEFit, dataset, CV(Ridge(lambda=0.1), 5, 42))` のまま、CV 1 タイプ追加だけで実現できる。MLJ.jl 流の "wrapping model" と同じパターン。

なお **λ-path**（resample なしで λ を変えるだけ）は `solve_coefficients` レベルで実装可能。つまり dispatch 位置が違うだけ:

| 種類 | dispatch レベル |
|---|---|
| λ-path（`RegularizationPath`） | `solve_coefficients` |
| k-fold CV / LOO-CV | `fit` |

### Bayesian 線形回帰の組み込み方

Bayesian は posterior（mean + covariance / samples）を返したいので、**result type を抽象化**する必要がある。GLM.jl の `fit(LinearModel, ...)` vs `fit(GeneralizedLinearModel, ...)` と同じパターンで:

```julia
abstract type AbstractSCEFit end

struct SCEFit <: AbstractSCEFit
    # 標準フィールド（dataset 参照、係数、メトリクス、残差）
end

struct BayesianSCEFit <: AbstractSCEFit
    # SCEFit のフィールド + posterior_mean / posterior_covariance（or samples）
end

fit(SCEFit, dataset, ::Bayesian) -> BayesianSCEFit  # 戻り値はサブタイプ
```

ユーザーコードは `fit(SCEFit, dataset, Bayesian(...))` のまま、返り値の型だけが豊かになる。アクセサ:
- `coef(fit)` → posterior mean
- `intercept(fit)` → posterior mean of bias
- `predict(SCEModel(fit), x)` → MAP / posterior mean prediction
- `posterior(fit)` / `samples(fit)` → 分布情報（Bayesian 固有、`SCEFit` には未定義）

`solve_coefficients(::Bayesian, X, y; bias_col)` の戻り値は j_values だけでは足りないので、NamedTuple 化（`(coefs=..., posterior=...)`）か、Bayesian の `fit` を solve_coefficients を経由せずに直接 override するか、いずれか。spec 化時に決める。

### まとめ

| 拡張対象 | 必要な追加 | API 影響 |
|---|---|---|
| λ-path | `solve_coefficients(::RegularizationPath, ...)` 1 メソッド | なし |
| k-fold CV | `fit(SCEFit, dataset, ::CV; ...)` 1 メソッド + `kfold` ユーティリティ | なし |
| Bayesian | `BayesianSCEFit <: AbstractSCEFit` + `fit(SCEFit, dataset, ::Bayesian)` + アクセサ | result type が hierarchy 化 |

本設計は CV と Bayesian に開かれているが、Bayesian を入れる段階で `AbstractSCEFit` を導入する必要があることだけ意識しておく。

## 既存 API との関係（コンビニエンス層として残すもの）

- `System(input_dict::Dict)` / `System(toml_file::String)`: pre-release のうちに削除。後方互換 alias は基本残さない（CLAUDE.md 方針: 数値・物理規約を変えない限り、API は理想形優先）。
- `SpinCluster`: 削除。新 API の `SCEDataset` + `SCEFit` で代替。
- `Optimizer` struct: 内部実装に格下げ。ユーザー向け result type は `SCEFit` / `SCEModel`。
- `write_xml`: `save(target, path)` のシノニムとして残すか deprecate するかは spec で決める（Julia の `JLD2.save` / `FileIO.save` との整合を考慮）。

## 設計判断ポイント（spec 化時に決める）

| 判断項目 | 候補 | 推奨と理由 |
|---|---|---|
| 物質型の名前 | `Crystal` / `Structure` を維持 | `Structure` 維持。既存コードへの影響最小化。`Crystal` は AtomsBase 等で使う語彙だが現状は導入しない。 |
| `fit` の引数順 | `fit(SCEFit, dataset, est)` / `fit(est, dataset)` | StatsAPI 準拠で type-first。GLM.jl の `fit(LinearModel, X, y)` と整合。 |
| `Symmetry` / `Cluster` の expose | 暗黙構築のみ / 明示構築も許可 | 明示構築も許可。default は暗黙、advanced 用途で `Symmetry(crystal; ...)` を直接構築可能に。 |
| Cross-validation の表現 | wrapper estimator `CV(...)` / 別関数 `crossvalidate(...)` | wrapper estimator を推奨。`fit` の API を変えず、estimator dispatch で吸収。 |
| AtomsBase.jl 統合 | 統合する / 当面なし | 当面なし。将来 `Crystal` から AtomsBase 型を生成する変換関数を足すだけで対応可能なので延期。 |
| 既存 `Optimizer` の扱い | export 維持 / 内部化 / 削除 | 内部化（または削除）。`SCEFit` / `SCEModel` がユーザー向け result type。 |
| `predict_torque` の置き場 | 独自関数 / `predict(...; what=:torque)` kwarg dispatch | 独自関数推奨（StatsAPI の `predict` は単一出力前提）。 |
| `SCEDataset` の field 公開度 | direct access / accessor 関数経由 | accessor 経由を推奨（`design_matrix_energy(dataset)` 等）。 |
| `weight` の表現 | kwarg `weight::Real` / 位置引数 | kwarg 推奨。デフォルト値（現状 `0.5`）と解釈の正規化方法（torque 行に `√(w/(1-w))` を掛ける等）も spec で明文化。 |
| StatsAPI.jl 依存追加 | 追加する / Magesty 内定義 | 軽量パッケージなので追加推奨。代替として `fit` / `coef` を Magesty 内定義で先行、後でアダプタ追加も可。 |

## 影響範囲（中規模 API リネーム + 拡張）

- 公開 export: `System`, `SpinCluster` 削除、`SCEProblem` / `SCEDataset` / `SCEFit` 追加（`SCEModel` 据置）。
- `build_sce_basis` / `build_sce_basis_from_xml` の戻り値型 → `SCEProblem`。
- `fit_sce_model(...)` → `fit(SCEFit, dataset, estimator; weight)` に置き換え。
- `write_xml(...)` → `save(model::SCEModel, path)` に置き換え（または `write_xml` をシノニムで残す）。
- ドキュメント全般（`SPEC.md`, `docs/src/`, README, `test/examples/`）。

## 移行戦略（spec 化したときの粒度）

1. `SCEProblem` 型と Julia kwargs constructor の追加（既存 `Structure`/`Symmetry`/`Cluster`/`BasisSet` を内部で構築する薄いラッパー）。
2. `SCEDataset` 型を新規導入（unweighted X_E/X_T を保持）。スライス / `vcat` API を実装。
3. estimator dispatch リファクタを完了（`docs/design-notes/estimator-dispatch.md` 参照、既に完了）。
4. `fit(SCEFit, dataset, estimator; weight)` の StatsAPI シグネチャを追加。
5. `coef` / `intercept` / `predict` / `residuals` / `metrics` を StatsAPI 準拠で実装。
6. `SCEModel(fit::SCEFit)` の lightweight 変換を実装、`save` / `load` を整理。
7. examples / test/examples を Julia スクリプトで書き直し。
8. `System` / `SpinCluster` / `fit_sce_model` を削除。

各段階は単独で merge 可能だが、最後の削除ステップは破壊的変更なので pre-release 中に確定させる。

## やらないこと（現時点）

- MLJ.jl / StatsModels.jl の正式統合（外部公開時に検討）。
- AtomsBase.jl 統合（将来オプション）。
- TOML の deprecation（コンビニエンス層として恒久的に残す）。
- input.toml フォーマット自体の改変（後方互換維持）。

## 進め方

`docs/specs/[YYMMDD]-sce-public-api/` を切って requirements.md / design.md / tasklist.md で合意してから着手する。estimator dispatch リファクタ（`docs/design-notes/estimator-dispatch.md`、完了済み）を前提として内部実装が簡潔になる。
