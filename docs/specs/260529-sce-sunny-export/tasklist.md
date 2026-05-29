# Tasklist: SCE → Sunny.jl LSWT エクスポータ

Status: in progress (2026-05-29) — M1/M3 done; M2 spike + M4/M5/M6 remaining

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

- [ ] `test/component/test_sunny_export.jl`（Sunny 非依存）。
- [ ] Sunny をテスト環境のみに追加 + ゲート設定（`MAGESTY_TEST_SUNNY`）。
- [ ] 往復検証テスト。既存レイアウトに合わせ fixture サブディレクトリ単位で配置
      （例 `test/integration/fept_tetragonal_2x2x2/`、`dimer/`）。各 fixture で
      `energy(sys) ≈ predict_energy - j0`。

### M6 — CLI・ドキュメント

- [ ] `cli/src/MagestyCLI.jl` に `@cast module Sunny`（`script`）＋ CLI テスト。
- [ ] `docs/src/sunny_export.md`、`docs/make.jl` ナビ、CLI ドキュメント、
      `docs/src/api.md`。

## Exit checklist

実装着地時に一度全項目を確認。該当しない項目は ~~取り消し線~~。

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] If results changed: regression or validation test added.（往復テスト）
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~（ホットパス非該当）
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
