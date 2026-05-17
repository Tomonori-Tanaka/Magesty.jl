# Pre-release Cleanup (v0.1.0)

**Status**: 未着手 (2026-05-17)

最初の正式リリース v0.1.0 に向けたパッケージ全体の点検結果と修正項目。examples パスが
壊れている (P0)、`Project.toml` の version/authors/compat 整理、`CITATION.cff` 追加、
Documenter 検査の強化、`predict_*` docstring 形式統一、`square_lattice` 未実装テストの
解消を一括で実施する。

## 背景

- 現在の `Project.toml` は `version = "0.1.0-DEV"` で、最新タグ `v0.0.2` から 60+ commits 経過。
- パッケージ全体の API・ドキュメント・テスト・CI を 3 エージェントで網羅調査し、リリースを
  阻む不具合と衛生問題を洗い出した。
- ユーザー回答: 「P0 修正 + メタデータ整理 + Doc/test 衛生 + version bump (0.1.0) + compat
  を minor まで緩める」をまとめて 1 PR 範囲に含める方針で合意。

## 検出された問題と修正方針

### P0. 壊れた examples のパス（即修正対象）

`FIXTURE` が `test/examples/` を指しているが実在するのは `test/integration/`。

- `examples/01_basic_flow.jl:10` — `"..", "test", "examples", "fept_tetragonal_2x2x2"` →
  `"..", "test", "integration", "fept_tetragonal_2x2x2"`
- `examples/03_save_load.jl:12` — 同上

`test/integration/fept_tetragonal_2x2x2/{input.toml,EMBSET.dat}` の存在は確認済み。
`examples/02_cif_input.jl` は問題なし。

### P1. Project.toml の整理

- **version** (行 3): `"0.1.0-DEV"` → `"0.1.0"`
- **authors** (行 4): `"T.Tanaka <tomonori.tanaka.academic@gmail.com>"` のまま据え置き
- **LICENSE 行 375**: `Copyright (c) 2024 TomonoriTanaka` を `T.Tanaka` に統一
- **compat の minor 化** (patch ピンを minor に緩和):
  - `AtomsBase = "0.5.2"` → `"0.5"`
  - `Combinat = "0.1.4"` → `"0.1"`
  - `EzXML = "1.2.1"` → `"1.2"`
  - `LegendrePolynomials = "0.4.5"` → `"0.4"`
  - `MultivariateStats = "0.10.3"` → `"0.10"`
  - `OffsetArrays = "1.17.0"` → `"1"`
  - `StaticArrays = "1.9.10"` → `"1.9"`
  - `Unitful = "1.28.0"` → `"1"`
  - `WignerD = "0.1.4"` → `"0.1"`
  - `JET = "0.11"` はテスト用 CI ピン目的で据え置き

`test/component/test_Version.jl` は `Project.toml` の version と `Magesty.VERSION` の文字列
一致を見るだけなので、bump で他テストへの波及はない。

### P2. CITATION.cff の追加（新規）

GitHub の "Cite this repository" が機能するように CFF 1.2.0 形式で。論文情報は README の
arXiv:2512.04458 / ORCID から取得。詳細フォーマット例:

```yaml
cff-version: 1.2.0
message: "If you use Magesty.jl, please cite the following work."
title: "Magesty.jl"
authors:
  - family-names: Tanaka
    given-names: Tomonori
    orcid: "https://orcid.org/0000-0001-7306-6770"
version: 0.1.0
date-released: 2026-05-17
license: MPL-2.0
repository-code: "https://github.com/Tomonori-Tanaka/Magesty.jl"
preferred-citation:
  type: article
  authors:
    - family-names: Tanaka
      given-names: Tomonori
  title: "General spin models from noncollinear spin DFT via spin-cluster expansion"
  year: 2025
  url: "https://arxiv.org/abs/2512.04458"
```

論文タイトル・著者リスト・arXiv ID は実装時に最終確認する。

### P3. Documenter 検査の強化

- `docs/make.jl:38` — `checkdocs = :none` → `checkdocs = :exports`
- 強化後、`export` シンボルすべてが `docs/src/api.md` 等で `@docs` ブロックを持つことを
  CI で検証できる。事前調査では網羅されていそうだが、ビルドが pass するかを実機で再確認。

### P4. predict_energy / predict_torque docstring の形式統一

- `src/Magesty.jl:529–540` (predict_energy)
- `src/Magesty.jl:559–571` (predict_torque)
- 現状は説明文のみ。CLAUDE.md の規約に従い `# Arguments` / `# Returns` を追加。同ファイル内の
  `r2_energy` / `rmse_torque` 等のスタイルに揃える。

### P5. 未実装テストの解消

- `test/integration/square_lattice/test.jl:85` — `# TODO: Add test for SCEModel` を解消。
- dimer / chain / febcc_2x2x2_pm 既存パターンを参考に、`SCEModel(fit)` 構築 +
  `predict_energy(model, dataset)` の往復チェックを追加する。`@test_skip` 化はしない。

### P6. オプション: tools/test/ の CI 統合

- `.github/workflows/CI.yml` に `make test-tools` step を追加。
- リリースをブロックしないが、補助スクリプトの腐敗を防ぐ。実施可否はユーザー判断。

## 検証手順

1. `julia --project=examples examples/01_basic_flow.jl` と `examples/03_save_load.jl` の完走確認。
2. `make test-all` — ユニット＋統合 pass（square_lattice / Version / Aqua 含む）。
3. `make test-jet` / `make test-aqua` — 警告 0。
4. `make test-tools` — tools スクリプトのテスト pass。
5. `julia --project=docs docs/make.jl` — `:exports` 強化後にビルド pass、未文書化 export 無し。
6. push 後 GitHub の "Cite this repository" ボタンが CITATION.cff を認識することを目視確認。

## 検出されなかった（対象外の）項目

事前調査で問題ないことを確認済み:

- `.jl` ソース中の日本語コメント・docstring
- `CLAUDE.md` / `docs/specs/` / `.claude/` への参照
- "spin cluster expansion" のハイフン抜け（src/Magesty.jl:5 で正しくハイフン有り）
- アメリカ英語規約違反（SpheriCart の `normalisation=:L2` のみ例外、正しく維持）
- TODO / FIXME / XXX / HACK コメント（square_lattice の 1 件のみ、P5 で対応）
- デバッグ用 `@show` / `println` 残骸
- `@test_broken` / `@test_skip` の不要残骸
- コメントアウトされたコードブロック
- `.DS_Store` の commit
- examples 02_cif_input.jl の API 整合性
- `docs/src/technical_notes.md` の数学規約と src 実装の整合
- 型アノテーション / `@views` / `@inbounds` のホットパス適用状況

## 進め方

- 本 design-note で合意 → spec 化するか直接 PR 化するかは別途判断。
- 実装時は機能単位で commit を分割（P0 / P1 / P2 / P3+P4 / P5 / 任意 P6）。
- `git commit` / `push` は user の明示指示後に `git-helper` agent 経由で実施。
- リリースタグ `v0.1.0` の作成は本作業 merge 後にユーザー側で実施想定。
