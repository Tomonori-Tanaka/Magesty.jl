# Magesty.jl — 仕様・アーキテクチャ

Julia で磁性材料の有効スピンモデルを構築するパッケージ。
ノンコリニアなスピン DFT 計算からスピンクラスター展開 (SCE) を用いて一般的なスピンモデルを構築する。

## 主要モジュール（`src/`）

| モジュール | ファイル | 役割 |
|-----------|---------|------|
| `Structures` | `src/Structures.jl` | 結晶構造の読み込み・スーパーセル生成 |
| `Symmetries` | `src/Symmetries.jl` | 空間群対称性操作（Spglib ラッパー） |
| `Clusters` | `src/Clusters.jl` | クラスター展開・距離行列計算 |
| `BasisSets` | `src/BasisSets.jl` | SALC 基底関数の生成（計算コスト大） |
| `Optimize` | `src/Optimize.jl` | スピン配置の最適化（OLS / Ridge） |
| `SpinConfigs` | `src/SpinConfigs.jl` | スピン配置の読み込み・管理 |

ユーティリティ（`src/utils/`）:

| ファイル | 役割 |
|---------|------|
| `SphericalHarmonicsTransforms.jl` | 実数球面調和関数の変換 |
| `MySphericalHarmonics.jl` | 球面調和関数の評価（unsafe 版含む） |
| `AngularMomentumCoupling.jl` | Wigner 係数・角運動量結合 |
| `ConfigParser.jl` | TOML 設定ファイルのパース |
| `RotationMatrix.jl` | 回転行列 |
| `xml_io.jl` | XML 形式での基底関数・結果の入出力 |
| `EnergyTorque.jl` | エネルギー・トルクの計算 |

共通データ構造（`src/common/`）: `SortedContainer`, `CountingContainer`

型定義（`src/types/`）: `AtomCells`, `AtomicIndices`, `Basis`

## 主要型

```
System
├── structure::Structure    # 結晶構造（ユニットセル・スーパーセル）
├── symmetry::Symmetry      # 対称操作・並進対称性
├── cluster::Cluster        # クラスター情報
└── basisset::BasisSet      # SALC 基底関数リスト

SpinCluster（System の拡張）
├── structure, symmetry, cluster, basisset（System と同じ）
└── optimize::Optimizer     # 最適化結果（SCE 係数・予測値）
```

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
examples/          利用例（FePt 等）
```

## 公開 API（`export` されているもの）

```julia
# 構築
System(input_dict)           # 基底関数まで構築
SpinCluster(input_dict)      # 最適化まで含めて構築
build_sce_basis(input_dict)  # 基底関数のみ構築
build_sce_basis_from_xml(input_dict, xml_file)  # XML から基底関数をロード

# 計算
calc_energy(sc, spin_config)   # エネルギー計算
calc_torque(sc, spin_config)   # トルク計算

# フィッティング
SCEModel, fit_sce_model, predict_energy
AbstractEstimator, OLS, Ridge

# パラメータ取得
get_j0(sc)        # 基準エネルギー J₀
get_jphi(sc)      # SCE 係数ベクトル
get_j0_jphi(sc)   # 両方

# 出力
write_xml(sc, filename)       # XML 出力（基底関数・係数）
write_energies(sc, filename)  # エネルギー一覧
write_torques(sc, filename)   # トルク一覧

# その他
read_embset(...)    # スピン配置の読み込み
install_tools()     # CLI ツールのインストール
```

## 主要な外部ライブラリ

| ライブラリ | 用途 |
|-----------|------|
| `Spglib` | 空間群解析 |
| `StaticArrays` | スタック割り当て配列（高速化） |
| `EzXML` | XML 入出力 |
| `WignerD`, `WignerSymbols`, `LegendrePolynomials` | 角運動量結合 |
| `LinearAlgebra`, `Statistics` | 線形代数・統計 |
| `MultivariateStats` | Ridge 回帰（`ridge`） |
