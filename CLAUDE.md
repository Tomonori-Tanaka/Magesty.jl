# CLAUDE.md

## プロジェクトの目的

ノンコリニアなスピン DFT 計算からスピンクラスター展開 (SCE) を用いて、一般的な有効スピンモデルを構築する Julia パッケージ。
スタイルのリファクタリングよりも、**数値的な正確さ・再現性・物理的な整合性**を優先する。

## 基本ルール

- 数値規約（符号・単位・規格化）を黙って変更しない。
- アルゴリズムを編集する前に、関連する方程式と現在の符号・単位規約を確認すること。
- 数値結果を変える可能性のある変更には、必ず以下を伴うこと：
  1. 結果が変わる理由の簡潔な説明
  2. リグレッションまたは検証テスト
  3. ユーザー向けの場合は `docs/` / `examples/` の更新
- `git commit` および `git push` は必ず事前に確認をとる（ローカルコミットも勝手に行わない）。

## 実装規約

- 隠れたグローバル状態を避ける。
- エクスポートされる API には明示的な型アノテーションと docstring を付ける。
- 型アノテーションを積極的に付ける（`::Float64`, `::Vector{Float64}`, `::SVector{3, Float64}` 等）。
- 公開 API の docstring は Julia 標準形式（`# Arguments` / `# Returns` / `# Examples`）で記述する。
- パフォーマンス改善は、変更前後でベンチマークを取り `DESIGN_NOTES.md` または `.claude/bench_log.md` に記録する。

詳細なコーディングスタイル（命名規則・ループ規約・引数順序など）は `STYLE_GUIDE.md` を参照。**コード編集時は必ず確認する**。

## 言語・用語規約

- **会話・説明文**: 日本語可。
- **ソース・コメント・docstring・コミットメッセージ・PR**: 英語。
- **コミットメッセージは [Conventional Commits](https://www.conventionalcommits.org/) に準拠する**。
  形式: `<type>(<scope>): <subject>`（scope は省略可）。
  - 使用する type: `feat` / `fix` / `docs` / `test` / `refactor` / `perf` / `chore` / `style`。
  - subject は命令形・小文字始まり・末尾ピリオドなし（例: `add SCE basis loader`）。
  - 破壊的変更は body に `BREAKING CHANGE: ...` を付ける。
- **アメリカ英語**で統一する（`normalize` / `behavior` / `color` / `center` 等）。
- 例外: **外部 API のリテラルは元の綴りを保つ**（例: SpheriCart の `normalisation=:L2` kwarg）。
- **`.jl` ソース中の日本語は PostToolUse hook (`.claude/hooks/no-japanese.sh`) で自動ブロック**。対象は `src/` `test/` `tools/` 配下（`tools/personal/` は除外）。
- **ソースコード中に Claude 内部ドキュメントを参照しない**: `.jl` のコメント・docstring 内で `CLAUDE.md` / `DESIGN_NOTES.md` / `docs/design-notes/` / `.claude/` / `docs/specs/` を**名指しで参照しない**。これらは Claude との協働用 scaffolding であり、公開ソースは独立して読めるべき。歴史的な経緯を残したい場合は内容をインラインに要約する。
  - 参照してよい外部ドキュメント: `docs/src/` (Documenter)、`SPEC.md`、`STYLE_GUIDE.md`、`https://Tomonori-Tanaka.github.io/Magesty.jl/` (technical notes).
  - 例外: コミットメッセージの `Refs:` 行で `docs/design-notes/refactor-sweep.md R7` のように参照するのは可（commit は workflow artifact）。

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

ベンチマーク用のターゲット（`make bench-sphericart` / `make bench-optimize-sphericart` 等）は Makefile を参照。

## 物理規約

これらを誤解すると無言でバグが埋まるので、アルゴリズムを触る前に確認すること。

- **スピン方向は単位ベクトル**: `spin_directions[:, atom]` は古典スピンの**方向**（`‖·‖ = 1`）。磁気モーメントの大きさは別概念。
- **スピン行列のレイアウト**: `3 × n_atoms`（行 = x,y,z、列 = 原子）。転置すると全計算が壊れる。
- **球面調和関数は実数（テッサー型）**: 内部表現は実 `Zₗₘ`（複素 `Yₗₘ` ではない）。`(l, m)` ペア数は `(l_max+1)²`。
  - ホットパスは `Zₗₘ_unsafe`（境界チェックなしの buffered 版）。意味は同じだが速い。
- **SALC・CG 係数の規約は本リポジトリで定義**: 変更前は必ず [Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/) を参照すること。
- **エネルギー単位**: SCE 係数 `Jφ` の単位は DFT 入力のエネルギー単位に従う（通常 eV）。`j0`（基準エネルギー）は SCE 係数とは別に保持される。

## パフォーマンス最適化方針

ホットパスは `Fitting.jl` の基底関数評価ループと `SALCBases.jl` の SALC 計算。

**StaticArrays の活用**:
- 3 次元ベクトルは `SVector{3, Float64}`（不変）または `MVector{3, Float64}`（スクラッチバッファ）を使う。
- ループ内で新しい `Vector` を生成しない。`MVector` を事前確保して再利用する。
- `spin_directions[:, atom]` の列スライスは `@views` でコピーを避け、`SVector` に変換してスタック上で処理する。

```julia
# Good
dir_svec = SVector{3, Float64}(spin_directions[:, atom])
buf = MVector{3, Float64}(0.0, 0.0, 0.0)

# Bad（ヒープ割り当てが発生）
dir_vec = spin_directions[:, atom]
buf = zeros(3)
```

**スレッド並列化**:
- スピン配置単位のループ（`num_spinconfigs`）が主要な `@threads` 対象。
- 対称操作ループ（`isym in 1:n_operations`）も `@threads` で並列化可能。

**境界チェック抑制**:
- インデックスの正しさが保証できる内側ループには `@inbounds` を付ける。
- 球面調和関数のホットパスには `Zₗₘ_unsafe`（境界チェックなし版）を使う。

**計算コストの大きい処理**:
- `SALCBasis` の SALC 計算は重い。再利用する場合は `write_xml` で保存し `build_sce_basis_from_xml` でロードする。

## 連動箇所（一方を変えたら全箇所を確認）

### 球面調和関数の規約
`TesseralHarmonics.jl` の `Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` の規格化と符号は `SphericalHarmonicsTransforms.jl` / SALC 構築（`SALCBases.jl`） / SpheriCart 比較テスト（`test/component_test/test_sphericart_agreement.jl`）と整合している必要がある。片方だけ変えると design matrix が静かに壊れる。

### SCE 係数の入出力
`write_xml` / `build_sce_basis_from_xml` のラウンドトリップは整数桁まで一致しなければならない。係数や基底順序のフォーマットを変える場合は両方を同期する。

### Fitting と SALCBasis の対応
`Fitting.jl` の design matrix 構築は `SALCBasis` 内の `(l, m, site)` 順序に依存している。基底の並び順を変えると係数の物理的解釈が変わる。

## 開発単位の管理

中規模以上の開発は **spec フォルダ** に集約する。スプリント横断の進捗トラッキングは作らない。

- **進行中の開発単位**: [`docs/specs/[YYMMDD]-[slug]/`](docs/specs/)（requirements.md / design.md / tasklist.md の 3 ファイル構成）。
- **横断的な設計メモ・調査結果・保留アイデア**: [`DESIGN_NOTES.md`](DESIGN_NOTES.md)（インデックス）+ [`docs/design-notes/`](docs/design-notes/)（トピック別本体）。運用ルールは `docs/design-notes/README.md`。
- **日々の細かい作業**: `TaskCreate`（in-session のみ。永続化したい軽い TODO は `docs/design-notes/` か spec へ）。
- 実装済みのベンチマーク履歴: `.claude/bench_log.md` と `git log`。

### spec フォルダの運用ルール

**中規模以上の開発を始める前に、必ず先に spec フォルダを作って合意を取る。**

判定基準（どれかに当てはまったら spec を作る）:
- 数日以上かかる。
- 設計判断が複数ある（API・型・規約に影響する選択）。
- 既存挙動が変わる中規模以上の変更（新機能追加 / 中規模リファクタリング）。
- 後から「なぜこう作った？」と訊かれそう。

spec を作らなくてよいもの:
- バグ修正（テスト追加で完結）。
- ドキュメント・コメント修正。
- 1 ファイル内の小規模 refactor。
- 既存テストで担保される動作の小修正。

手順（Claude 側で実施）:
1. `docs/specs/[YYMMDD]-[slug]/` を作る（`YYMMDD` = `date +%y%m%d`、`slug` は英語の kebab-case）。
2. 同フォルダに以下 3 ファイルを置き、user と相談しながら埋める:
   - `requirements.md` — 目的・不変条件・スコープ・完了基準。
   - `design.md` — モジュール構成・API・型・規約・連動箇所への影響。
   - `tasklist.md` — マイルストーン（粗い粒度。日々の細かい作業は TaskCreate）。
3. spec の合意ができてから実装に着手する。
4. 完了後もフォルダは残す（削除しない、履歴として参照）。

## Claude の作業方針

### 確認不要（自由にやってよい）

- バグ修正（最小限の変更 + テスト追加）。
- テストの追加・修正。
- ドキュメントの誤記修正。
- `DESIGN_NOTES.md` インデックス / `docs/design-notes/` 配下への気付きの追記。

### サブエージェントの活用

- 実装後は `test-runner` エージェントでテストを実行・診断する（`.claude/agents/test-runner.md`）。
- コミット前は `code-reviewer` エージェントで変更差分をレビューする（`.claude/agents/code-reviewer.md`）。
- 性能調査時は `profiler` エージェントを使う（`.claude/agents/profiler.md`）。
- user から commit / push の明示指示を受けた後は `git-helper` エージェントを呼ぶ（`.claude/agents/git-helper.md`）。Conventional Commits 起草・ASCII 検査・`Write` + `git commit -F file` 経由のコミット適用・BREAKING CHANGE 自動検出までを一括で行う。heredoc 事故を避けるため、メインエージェントが直接 `git commit -m` を実行しない方針。

### 提案してから実装する

- アルゴリズムの変更（数値結果が変わる可能性があるとき）。
- リファクタリング（モジュール境界をまたぐもの）。
- パフォーマンス改善（先にベンチマーク結果を示す）。

### 実装せず、必ず確認する

- 物理規約の変更（符号・単位・規格化・SALC や CG 係数の規約）。
- 新しい外部依存の追加。
- 公開 API（`export` されているもの）のシグネチャ変更。
- 球面調和関数の規格化変更（`Zₗₘ` 系列）。
- XML 入出力フォーマットの変更。
- `git commit` / `git push` 全般（ローカルコミットを含めて、明示的な指示なしに実行しない）。

## 参照

作業前に、必要に応じて以下を参照する。

- `STYLE_GUIDE.md` — コーディングスタイルの詳細ルール（**コード編集時は必ず確認する**）。
- `SPEC.md` — アーキテクチャ・ディレクトリ構成・主要型・公開 API。
- `docs/specs/` — 進行中・完了済みの spec（中規模以上の開発単位）。
- `DESIGN_NOTES.md` — 設計メモ・調査結果・保留アイデアのインデックス。本体は `docs/design-notes/` 配下。
- [Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/) — SALC・CG 係数等の数学的規約。
