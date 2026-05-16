# Magesty.jl — 仕様・アーキテクチャ

Julia で磁性材料の有効スピンモデルを構築するパッケージ。
ノンコリニアなスピン DFT 計算からスピンクラスター展開 (SCE) を用いて一般的なスピンモデルを構築する。

## 主要モジュール（`src/`）

すべてのモジュールは `src/` 直下にフラットに配置している。`include`
順序は `src/Magesty.jl` を参照（依存関係の順）。

### SCE パイプライン

| モジュール | ファイル | 役割 |
|-----------|---------|------|
| `Structures` | `src/Structures.jl` | 結晶構造の読み込み・スーパーセル生成 |
| `Symmetries` | `src/Symmetries.jl` | 空間群対称性操作（Spglib ラッパー） |
| `Clusters` | `src/Clusters.jl` | クラスター展開・距離行列計算 |
| `SALCBases` | `src/SALCBases.jl` | SALC 基底関数の生成（計算コスト大） |
| `Fitting` | `src/Fitting.jl` | デザイン行列構築・回帰（OLS / Ridge） |
| `SpinConfigs` | `src/SpinConfigs.jl` | スピン配置の読み込み・管理 |

### ユーティリティ・データ型

| モジュール / ファイル | 役割 |
|---------|------|
| `SphericalHarmonicsTransforms.jl` | 実数球面調和関数の変換 |
| `TesseralHarmonics.jl` | テッセラル球面調和関数 `Zₗₘ` の評価（unsafe 版含む） |
| `AngularMomentumCoupling.jl` | Wigner 係数・角運動量結合 |
| `CoupledBases.jl` | 角運動量結合で組んだ coupled basis の中間表現 |
| `InputSpecs.jl` | 入力 typed value (`SystemSpec` / `InteractionSpec` / `SymmetryOptions`) と TOML / Dict パーサ (`parse_toml_inputs`、ワイルドカード species 展開を含む) |
| `RotationMatrix.jl` | 回転行列 |
| `XMLIO.jl` | XML 形式での基底関数・モデルの入出力 |
| `AtomsBaseAdapter.jl` | AtomsBase / Unitful の境界アダプタ (`system_to_specs` / `kwargs_to_specs`) |
| `AtomCells.jl` | 原子サイト + ユニットセル情報の軽量データ型 |
| `SortedCounters.jl` | ソート済みキー反復をサポートする内部カウンタ |
| `Version.jl` | Project.toml と整合するバージョン情報 |

## 主要型（新 API）

```
# 入力の typed value gate（src/InputSpecs.jl）。4 つの入力経路
# (TOML file / Dict / AtomsBase / kwargs) はすべてこの 3-tuple を組み立て、
# それを SCEBasis 内部コンストラクタに渡す。

SystemSpec
├── name::String
├── num_atoms::Int
├── kd_name::Vector{String}           # 種類名（順序が kd_int_list の番号付け）
├── kd_int_list::Vector{Int}          # 原子ごとの種類インデックス
├── lattice_vectors::Matrix{Float64}  # 3×3 (Å)
├── x_fractional::Matrix{Float64}     # 3×num_atoms
└── is_periodic::Vector{Bool}         # 長さ 3

InteractionSpec
├── nbody::Int
├── body1_lmax::Vector{Int}                  # 種類数の長さ
├── bodyn_lsum::OffsetArray{Int, 1}          # axis 2:nbody
└── bodyn_cutoff::OffsetArray{Float64, 3}    # (2:nbody, 1:nkd, 1:nkd), 対称

SymmetryOptions
├── tolerance_sym::Float64
└── isotropy::Bool

# 下流の SCE パイプライン

SCEBasis
├── structure::Structure    # 結晶構造（ユニットセル・スーパーセル）
├── symmetry::Symmetry      # 対称操作・並進対称性
├── salcbasis::SALCBasis    # SALC 基底関数リスト
└── isotropy::Bool          # 等方性制限の構築フラグ

SCEDataset
├── basis::SCEBasis
├── spinconfigs::Vector{SpinConfig}
├── X_E::Matrix{Float64}    # エネルギーのデザイン行列（unweighted、列 1 はバイアス）
├── X_T::Matrix{Float64}    # トルクのデザイン行列（unweighted、バイアス列なし）
├── y_E::Vector{Float64}    # 観測エネルギー
└── y_T::Vector{Float64}    # 観測トルク（flatten）

SCEFit <: StatsAPI.RegressionModel
├── dataset::SCEDataset
├── j0::Float64             # 推定された基準エネルギー
├── jphi::Vector{Float64}   # 推定された SCE 係数
├── estimator::AbstractEstimator
├── torque_weight::Float64
├── residuals::Vector{Float64}
└── metrics::Dict{Symbol, Any}  # in-sample RMSE / R²

SCEModel
├── basis::SCEBasis
├── j0::Float64
└── jphi::Vector{Float64}
```

`SCEModel(f::SCEFit)` で `SCEFit` から軽量な予測子に変換できる。

## ディレクトリ構成

```
src/               パッケージ本体
test/
  component_test/  ユニットテスト（モジュール単位）
  examples/        統合テスト（実際の計算例）
tools/             パッケージ本体とは独立したスクリプト群
  vasp/            VASP 入出力変換ツール
  test/            tools のテスト（make test-tools で実行）
  personal/        個人用スクリプト（品質保証対象外）
docs/              Documenter.jl ドキュメント
examples/          利用例（basic_flow / CIF input / save-load）
```

## 公開 API（`export` されているもの）

### 新 API

```julia
# 構築
SCEBasis(toml_path)                # TOML から SCE 基底を構築（SALC 含む）
SCEBasis(input_dict)               # parsed TOML dict から
SCEBasis(system; interaction, ...) # AtomsBase.AbstractSystem から
SCEBasis(; lattice, kd, kd_list, positions, periodicity, interaction, ...)

SCEDataset(basis, spinconfigs)     # SpinConfig ベクトルから
SCEDataset(basis, embset_path)     # EMBSET ファイルから
SCEDataset(system, spinconfigs; interaction, ...)  # AtomsBase shortcut
SCEDataset(toml_path, spinconfigs) # TOML shortcut

# フィッティング
fit(SCEFit, dataset, estimator; torque_weight = 0.5) -> SCEFit
SCEModel(f::SCEFit) -> SCEModel    # 予測専用の軽量変換

# 予測（model / fit / dataset / SpinConfig / Matrix を受け付ける）
predict_energy(model_or_fit, data)
predict_torque(model_or_fit, data)

# 評価（StatsAPI 規約 + Magesty ネイティブ）
coef(f), intercept(f), nobs(f), dof(f)
r2_energy(f), r2_torque(f)
rmse_energy(f), rmse_torque(f)
rss_energy(f), rss_torque(f)
residuals_energy(f), residuals_torque(f)

# 永続化（XML）
save(basis_or_model, path)         # path は .xml
load(SCEBasis, path)               # SCEBasis を読み出し（model XML でも可）
load(SCEModel, path)               # SCEModel を読み出し

# 推定器
AbstractEstimator, OLS, Ridge

# データ読み込み
read_embset(path)                  # EMBSET.dat -> Vector{SpinConfig}
```

## 主要な外部ライブラリ

| ライブラリ | 用途 |
|-----------|------|
| `AtomsBase` / `Unitful` | AtomsBase 入力経路 |
| `StatsAPI` | `fit` / `coef` / `nobs` / `dof` の規約 |
| `Spglib` | 空間群解析 |
| `StaticArrays` | スタック割り当て配列（高速化） |
| `EzXML` | XML 入出力 |
| `WignerD`, `WignerSymbols`, `LegendrePolynomials` | 角運動量結合 |
| `LinearAlgebra`, `Statistics` | 線形代数・統計 |
| `MultivariateStats` | Ridge 回帰（`ridge`） |
