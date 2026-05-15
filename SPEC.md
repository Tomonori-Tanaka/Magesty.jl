# Magesty.jl — 仕様・アーキテクチャ

Julia で磁性材料の有効スピンモデルを構築するパッケージ。
ノンコリニアなスピン DFT 計算からスピンクラスター展開 (SCE) を用いて一般的なスピンモデルを構築する。

## 主要モジュール（`src/`）

| モジュール | ファイル | 役割 |
|-----------|---------|------|
| `Structures` | `src/Structures.jl` | 結晶構造の読み込み・スーパーセル生成 |
| `Symmetries` | `src/Symmetries.jl` | 空間群対称性操作（Spglib ラッパー） |
| `Clusters` | `src/Clusters.jl` | クラスター展開・距離行列計算 |
| `SALCBases` | `src/SALCBases.jl` | SALC 基底関数の生成（計算コスト大） |
| `Optimize` | `src/Optimize.jl` | デザイン行列構築・回帰（OLS / Ridge） |
| `SpinConfigs` | `src/SpinConfigs.jl` | スピン配置の読み込み・管理 |

ユーティリティ（`src/utils/`）:

| ファイル | 役割 |
|---------|------|
| `SphericalHarmonicsTransforms.jl` | 実数球面調和関数の変換 |
| `MySphericalHarmonics.jl` | 球面調和関数の評価（unsafe 版含む） |
| `AngularMomentumCoupling.jl` | Wigner 係数・角運動量結合 |
| `ConfigParser.jl` | TOML 設定ファイルのパース |
| `RotationMatrix.jl` | 回転行列 |
| `xml_io.jl` | XML 形式での基底関数・モデルの入出力 |
| `EnergyTorque.jl` | エネルギー・トルクの計算 |
| `atomsbase_adapter.jl` | AtomsBase / Unitful の境界アダプタ |

共通データ構造（`src/common/`）: `SortedCounter`, `Version`

型定義（`src/types/`）: `AtomCells`, `Basis`

## 主要型（新 API）

```
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
  develop_tmp/     開発中・実験的テスト（CI 対象外）
  helpers/         テスト補助ユーティリティ
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

### レガシー API（将来削除予定）

```julia
System(input_dict)
SpinCluster(input_dict)
build_sce_basis(input_dict)
build_sce_basis_from_xml(input_dict, xml_file)
fit_sce_model(system, configs, est, weight)
write_xml(sc, path)
Magesty.calc_energy(sc, spin_config)
Magesty.calc_torque(sc, spin_config)
Magesty.get_j0(sc), Magesty.get_jphi(sc), Magesty.get_j0_jphi(sc)
Magesty.write_energies(sc, path), Magesty.write_torques(sc, path)
```

新 API への移行表は `docs/specs/260514-sce-public-api/design.md` の
"Migration map" を参照。

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
