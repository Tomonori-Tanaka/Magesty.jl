## コーディング規則

- **説明・コミットメッセージ**: 日本語
- **コードコメント**: 英語
- 型アノテーションを積極的に付ける（`::Float64`, `::Vector{Float64}` 等）
- 公開 API の docstring は Julia 標準形式（`# Arguments`, `# Returns`, `# Examples`）で記述

## テスト

修正後は必ず Makefile 経由でテストを実行する。

| コマンド | 対象 | 用途 |
|---------|------|------|
| `make test-unit` | `test/component_test/` | モジュール単位のユニットテスト |
| `make test-integration` | `test/examples/` | 実際の計算例を用いた統合テスト |
| `make test-all` | 上記両方 | 通常はこれを使う |
| `make test-tools` | `tools/test/` | tools スクリプトのテスト |
| `make test-jet` | — | JET.jl 静的型解析 |
| `make test-aqua` | — | Aqua.jl パッケージ品質チェック |

`test/develop_tmp/` は開発中の実験的テスト（CI 対象外）。

## パフォーマンス最適化方針

ホットパスは `Optimize.jl` の基底関数評価ループと `BasisSets.jl` の SALC 計算。

**StaticArrays の活用**:
- 3 次元ベクトルは `SVector{3, Float64}`（不変）または `MVector{3, Float64}`（スクラッチバッファ）を使う
- ループ内で新しい `Vector` を生成しない。`MVector` を事前確保して再利用する
- `spin_directions[:, atom]` の列スライスは `@views` でコピーを避け、`SVector` に変換してスタック上で処理する

```julia
# Good
dir_svec = SVector{3, Float64}(spin_directions[:, atom])
buf = MVector{3, Float64}(0.0, 0.0, 0.0)

# Bad（ヒープ割り当てが発生）
dir_vec = spin_directions[:, atom]
buf = zeros(3)
```

**スレッド並列化**:
- スピン配置単位のループ（`num_spinconfigs`）が主要な `@threads` 対象
- 対称操作ループ（`isym in 1:n_operations`）も `@threads` で並列化可能

**境界チェック抑制**:
- インデックスの正しさが保証できる内側ループには `@inbounds` を付ける
- 球面調和関数のホットパスには `Zₗₘ_unsafe`（境界チェックなし版）を使う

**計算コストの大きい処理**:
- `BasisSet` の SALC 計算は重い。再利用する場合は `write_xml` で保存し `build_sce_basis_from_xml` でロードする

## タスク管理

実装タスクと進捗は `TASK.md` で管理する。

## 改善案メモ

作業中に気付いた改善案・設計上の検討事項は `DESIGN_NOTES.md`（リポジトリ直下）に追記する。即座に対応しないアイデアもここに残しておく。

## 参照

作業前に、必要に応じて以下を参照してください。

- `STYLE_GUIDE.md` — コーディングスタイルの詳細ルール（**コード編集時は必ず確認する**）
- `SPEC.md` — アーキテクチャ・ディレクトリ構成・主要型の詳細
- `TASK.md` — 実装タスクと進捗
