# Tasklist: sce_to_sunny に物理スピン（実効スピン長）を渡す

Status: implemented (2026-06-18) — tests green (unit / sunny / jet / aqua) and
Tier 2 review panel resolved; pending commit on `feat/sunny-physical-spin`.

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — spin/g 解決とモード選択（コア・Sunny 非依存）

- [ ] `_sunny_resolve_spins(model, spin, g)`: スカラ／Dict、元素照合、非磁性
      placeholder、磁性種欠落エラー、`spin=nothing` の `ArgumentError`。
- [ ] `_sunny_onsite_factor(s, mode)`: `:dipole`→`2/(s(2s-1))`、
      `:dipole_uncorrected`→`1/s²`、`s=1/2`+`:dipole` ガード。
- [ ] `_sunny_select_mode(spins, mode)`: `:auto` の半整数判定。
- [ ] `_SUNNY_SPIN_S` 定数の撤去。

### M2 — emission での再スケール

- [ ] primitive ルート: ボンド `/(s_i s_j)`、onsite 係数化、`Moment` per-site。
- [ ] explicit ルート: 同上（スーパーセル原子単位）。
- [ ] ヘッダコメント更新（物理スピン・モード・`energy` 不変性）。

### M3 — テスト

- [ ] 単体（`test/component/test_sunny_export.jl`）: 解決・再スケール・onsite 係数・
      モード自動・エラー系・構文妥当性（既存 `placement` エラーテストは残す）。
- [ ] 往復（`test/sunny/runtests.jl`, `make test-sunny`）: `energy ≈ predict_energy - j0`
      を `spin × mode` 横断。
- [ ] 分散スケール回帰（同上）: `s=1` vs `s=5/2` のバンド幅比 ≈ 2.5。

### M4 — CLI・ドキュメント

- [ ] `magesty sunny script` に `--spin`/`--g`/`--mode`。
- [ ] `docs/src/tools.md`（`magesty sunny script` 節, narrative）＋ `api.md`：
      `S_eff = m/(gμ_B)`、`s(2s-1)/2`、モード選択、遍歴磁性の注意、文献。
- [ ] `SPEC.md` シグネチャ反映（exit checklist の「public API 変更」項と対応）。
- [ ] `CHANGELOG.md` `[Unreleased]` に `BREAKING CHANGE`。

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes. (`make test-unit` 22817 passed)
- [x] `make test-sunny` passes (Sunny-gated round-trip + dispersion-scale
      regression; not part of `test-all`). (27 passed)
- [x] `make test-aqua` / `make test-jet` clean (no new warnings). (aqua 10 / jet 1)
- [x] If results changed: regression or validation test added.
      （マグノン分散変化 → `1/s` スケール回帰で固定）
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
      （narrative `docs/src/tools.md` も更新）
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~（SunnyExport は hot path でない）
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [x] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~（不変の見込み）
- [x] `CHANGELOG.md` `[Unreleased]` updated（`BREAKING CHANGE`）.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.（commit 時に complete へ）
- [ ] Implementation commit hash appended below.
