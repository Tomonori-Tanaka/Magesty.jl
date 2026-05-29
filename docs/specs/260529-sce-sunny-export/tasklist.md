# Tasklist: SCE → Sunny.jl LSWT エクスポータ

Status: complete (2026-05-29) — pending merge of feature/sce-sunny-export

このファイルは commit 粒度の粗いマイルストーンを保持する。日々の作業は
セッション内の `TaskCreate` で追う。

## Milestones

### M1 — spec 合意

- [x] requirements / design / tasklist をユーザと合意。

### M2 — 早期 Sunny スパイク（完了）

- [x] エネルギー規約を実測確定（exchange は M 直接、onsite は `2/(s(2s-1))·ΣA SαSβ`、
      s=1→2A）。dimer/febcc/fept/fege で `:explicit` が機械精度一致。
- [x] reshape ショートカット不可を確認（Sunny は大スーパーセルで対称性解析を無効化）。
- [x] Magesty の spglib プリミティブから小セルを構築 → Sunny が完全対称性同定
      （febcc→Im-3m, fept→P4/mmm, fege→P2₁3）。
- [x] プリミティブ展開の正確な条件を確定: `multiplicity × #clusters == ntran`
      （cutoff ≤ L/2）。正しい実装（素のボンド結合・最小像オフセット・P1 強制）で
      **dimer・chain（cutoff=L/2 多セル）を機械精度一致**で実証。
- [x] placement 既定を `:primitive`（クリーン時）/ `:explicit`（フォールバック）に確定。

### M3 — 分解コア（`_sunny_decompose` / `_sunny_build_terms`）

注: 分解は placement に依存しないため M2 より先に実装・検証した。

- [x] `src/SunnyExport.jl` と IR（`_SunnyTerms`: pairs / onsites /
      const_offset / skipped）。
- [x] SALC 走査・分類・`folded_tensor` からの 3×3 構築・正準ボンド集約（転置含む）。
- [x] 単サイト `l=2` → 単イオン `A`（対称トレースレス）。
- [x] スキップ判定と警告収集、定数 SALC のオフセット畳み込み。
- [x] component テスト（規約突合 + dimer/fept/fege 再構成を機械精度で検証、
      skip ロジック）。`make test-unit` 緑。

### M4 — スクリプト生成・公開（完了）

- [x] プリミティブ構築（`_sunny_build_primitive`: spglib プリミティブ抽出＋素の
      ボンド結合＋最小像オフセット＋`mult×#clusters==ntran` クリーン判定）。
- [x] emitter（primitive/explicit 両ルート、フルパイプライン＋磁気構造ユーザ編集
      セクション＋スキップ警告）。
- [x] 公開 `sce_to_sunny(model; output, placement=:auto, symprec)`（自動ルート選択
      ＋警告フォールバック、docstring 付き）、`src/Magesty.jl` で export。
- [x] CI 単体テスト（クリーン判定・ルート選択・スクリプト parse・ファイル出力・
      引数検証）。`make test-unit` 緑。
- [x] **emitter end-to-end 検証**（生成スクリプトを eval して Sunny `energy(sys)` と
      `predict_energy-j0` を突合）: dimer/chain(primitive) 機械精度、
      febcc/fept(explicit, 単イオン含む) 〜1e-13。

### M5 — テスト

- [x] `test/component/test_sunny_export.jl`（Sunny 非依存、CI: 277 件）。
- [x] Sunny を専用サブ環境 `test/sunny/`（独自 Project.toml）に隔離。`make test-all`
      には含めず、`make test-sunny` ターゲットで実行（Sunny は重いため）。
- [x] 往復検証 `test/sunny/runtests.jl`：生成スクリプトを eval し `energy(sys)` と
      `predict_energy - j0` を突合。primitive（dimer/chain）・explicit（febcc/fept,
      単イオン含む）両ルート。`make test-sunny` 緑（15 件）。
- [ ] （任意・将来）3D のクリーン fixture（cutoff<L/2、複数ボンド）で primitive を
      さらに強化。chain は擬1D・単一ボンドのため。

### M6 — CLI・ドキュメント（完了）

- [x] `cli/src/MagestyCLI.jl` に `@cast module Sunny`（`script` サブコマンド）。
      `cli/test/runtests.jl` にテスト追加、`make test-cli` 緑。
- [x] ドキュメント：`docs/src/api.md`（`@docs sce_to_sunny` ＋ Sunny export 節）、
      `docs/src/tools.md`（`magesty sunny script` の解説＋変換表＋セル方針）、
      `SPEC.md`（モジュール表＋公開 API）、`CHANGELOG.md` [Unreleased] Added。
      （既存の `tools.md` に統合したため新規 `sunny_export.md` は作らず、make.jl ナビ
      変更も不要。）
- [x] `.claude/agents/test-runner.md` に `make test-sunny` を追記。
- [x] `make test-aqua`（10/10）・`make test-jet`（no issues）通過。

## Exit checklist

実装着地時に一度全項目を確認。該当しない項目は ~~取り消し線~~。

- [x] `make test-all` passes（unit 277 含む）。
- [x] `make test-aqua`（10/10）/ `make test-jet`（no issues）clean。
      ＋ `make test-sunny`（15、往復一致）。
- [x] If results changed: regression or validation test added（往復テスト
      `test/sunny/` ＋ component 再構成テスト）。
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated。
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~（ホットパス非該当）
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved. Numerical: 0 blocker/major（変換式・スケーリング・
      符号・単位・s=1 係数・DM フィルタ全て clean 確認）。適用: 実行時スキップ
      `@warn`、未使用 `symprec` 削除、CLI placement 検証、死コード `cellof` 削除、
      `:explicit` 時の二重 SALC 走査回避、s/係数の定数化、最小像 floor 丸め、型注釈・
      行長・コメント整備。
- [x] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated（test-runner に test-sunny 追記）。
- [x] `CHANGELOG.md` `[Unreleased]` updated。
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync。
- [x] Implementation commits appended below.

## Implementation commits (branch feature/sce-sunny-export)

- `30db236` chore: allow Japanese in per-slug spec working files
- `1aa145c` feat(tools): add SCE-to-Sunny exporter decomposition core
- `de6cb6d` feat(tools): emit Sunny.jl LSWT scripts via sce_to_sunny
- `7c86e2d` feat(cli): add 'magesty sunny script' and Sunny round-trip test
- `330712d` refactor(tools): apply Tier 2 review fixes to sce_to_sunny
