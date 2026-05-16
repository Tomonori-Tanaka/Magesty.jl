# Requirements: `src/` layout refactor for publication readiness

## Background

`v0.1.0-DEV` の公開（GitHub / General registry 登録を想定）に向けて、
`src/` 配下のファイル名・モジュール名・ディレクトリ構成を「初見の読み手」
にとって違和感なく追えるように整える。実装挙動は変えない — 名前と配置のみ。

## Goals

1. ファイル名とその中で `module X` 宣言する `X` を一致させる
   （Julia の標準的な慣習）。
2. `My` プレフィックスのように公開パッケージとして気が引ける名前を改める。
3. `common/` `types/` `utils/` の境界が曖昧でファイル数も偏っている
   現状を解消する。
4. `src/Magesty.jl` の `include` ブロックを上から読んだときに依存順序が
   理解しやすくなっていること。

## Non-goals

- 公開 API（`export` されているシンボル）のシグネチャ変更。
  改称対象は **内部モジュール名のみ**。
- 物理規約（符号・単位・規格化）の変更。
- パフォーマンス変更。
- テストの追加・修正以外でのテスト構造変更。
- `tools/` の構造変更（中身のリファクタは別件）。

## Invariants

R1. **公開 API は完全に不変**:
    `SCEBasis`, `SCEDataset`, `SCEFit`, `SCEModel`, `OLS`, `Ridge`,
    `AbstractEstimator`, `fit`, `coef`, `intercept`, `nobs`, `dof`,
    `predict_energy`, `predict_torque`, `r2_*`, `rmse_*`, `rss_*`,
    `residuals_*`, `save`, `load`, `read_embset`, `VERSION`,
    `install_tools` — すべて同じ呼び出しで動く。

R2. **テストは現状通り全て緑**:
    `make test-unit` / `make test-integration` / `make test-aqua` /
    `make test-jet` / `make test-tools` がリファクタ前と同じ結果。

R3. **XML スキーマと EMBSET フォーマットは無変更**:
    既存の `.xml` / `EMBSET.dat` がそのまま読み込める。

R4. **docs ビルドが通る**: `make build` がエラーなし。
    `@docs` ブロックで参照しているモジュール／関数のパスが追従する
    （e.g. `MySphericalHarmonics.Zₗₘ` → 改称後の名前）。

R5. **`SPEC.md` がリファクタ後のレイアウトを正しく反映**。

## Scope

### In-scope

- `src/` 配下のファイル名・ディレクトリ構成・内部モジュール名の変更
- `src/Magesty.jl` 内の `include` パスと `using .X` の調整
- `SPEC.md` のディレクトリ図と表の同期
- `docs/src/api.md` の `@docs` ブロック内モジュール名の同期
- `CLAUDE.md` 内で `MySphericalHarmonics` 等を名指ししている箇所の同期
- `test/` の `using` 文の更新（モジュール名が変わった場合）
- **`tools/` 配下の追随修正のみ**: `src/` の改名・移動に伴う
  `include`／`using`／パス参照の最小限の同期（構造の変更や中身の
  リファクタは含めない）

### Out-of-scope

- `tools/` の構造変更・スクリプトのリファクタ（別件）
- `docs/specs/`（過去 spec の本文は履歴として保持）
- `examples/`（公開 API のみ使うので無変更）

## Completion criteria

C1. `src/` 配下の `.jl` ファイル名と `module` 宣言が一致している。
C2. `My` プレフィックス等、再考対象として挙がった名前が全て解消されている。
C3. `make test-all` / `make test-aqua` / `make test-jet` /
    `make test-tools` / `make build` すべて緑。
C4. 公開 API テスト（`test/component_test/`）の `using Magesty` 起点の
    呼び出しはコード変更なしで動く。
C5. `SPEC.md` と `docs/src/api.md` がリファクタ後の構成を反映している。
C6. `tools/` 配下に `src/` の旧パス／旧モジュール名を指す参照が残って
    いない（`grep` 確認）。
C7. PR description に「公開 API 不変」と「改称マッピング表」が含まれる。
