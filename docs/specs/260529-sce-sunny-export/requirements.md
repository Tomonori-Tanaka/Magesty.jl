# Requirements: SCE → Sunny.jl LSWT エクスポータ

Status: draft (2026-05-29)

## Goal

フィット済みの `SCEModel` から、Sunny.jl で線形スピン波理論（LSWT）の
スピン波分散を計算する**実行可能な Julia スクリプト**を生成する新ツール
`sce_to_sunny` を追加する。DFT → SCE フィット → マグノン分散までを、係数の
手動翻訳なしにつなぐ。

## Background

- SCE モデルのエネルギーは `E = j0 + Σ_ν jphi[ν]·Φ_ν({ê_i})`。`ê_i` は単位
  スピン方向（`3 × num_atoms`）。
- 最低次の SCE 係数から従来型スピンモデル（Heisenberg `J` / DM ベクトル `D` /
  異方的対称テンソル `Γ`）への変換式は `docs/src/technical_notes.md` に導出・
  検証済み。本ツールはこれを実装に落とす。
- Sunny.jl は LSWT・動的構造因子の成熟した実装を持つため、Magesty 側は係数を
  与えるスクリプトを書き出すだけでよい。
- 既存の `vasp_to_extxyz` と**同じパターン**（`VaspConvert.jl` のように
  モジュールにせず Magesty 名前空間へ直接 include、内部ヘルパは `_sunny_*`
  プレフィックス、純関数で文字列を返し任意でファイル書き出し）に揃える。
  テキスト生成のみなので **Magesty コアに Sunny.jl を依存追加しない**。

## Scope

Includes:

- 新サブモジュール `src/SunnyExport.jl` と公開関数
  `sce_to_sunny(model::SCEModel; ...)::String`。
- 対応する相互作用（v1）:
  - 2体・`l1=l2=1` のペア項 → 3×3 交換テンソル（`set_exchange!`）。
  - 単サイト `l=2` 項 → 単イオン異方性（`set_onsite_coupling!`）。
- 上記以外（高次 `l` ペア、3体以上）は**検出してスキップ**し、何を落としたか
  を `@warn` とスクリプト内コメントで明示。
- フルパイプラインのスクリプト生成: Crystal → System → 相互作用 →
  `reshape_supercell`（プリミティブ化）→ `minimize_energy!` →
  `SpinWaveTheory` → `q_space_path` → `dispersion` → プロット。
- CLI サブコマンド `magesty sunny script`。
- Sunny による往復エネルギー検証テスト（テスト環境のみ、CI ゲート付き）。

Excludes:

- 物理的スピン長 `S` によるリスケール（既定は `s=1` 規約単位）。将来フックのみ
  用意。
- 高次相互作用（biquadratic、3体）の Sunny への明示エクスポート。
- 高対称点の自動決定（q-path はユーザ編集前提のプレースホルダ）。

## Invariants

- `sce_to_sunny` が表す相互作用での Sunny エネルギーは、定数項 `j0` を除き
  `Magesty.predict_energy` と一致する（往復テストで担保）。これが最重要不変条件。
- スピン方向は単位ベクトル・`3 × num_atoms` レイアウトのまま。
- 球面調和の規約・SALC ↔ Fitting の整合は変更しない（読み取りのみ）。
- Magesty コアの依存（`Project.toml` `[deps]`）に Sunny を追加しない。
- 既存の言語ポリシー: 本 spec 実体ファイルは日本語可、生成スクリプトやソースは
  英語。

## Completion criteria

- [ ] `make test-unit` が通る（Sunny 非依存の分解・スキップ・文字列検証）。
- [ ] Sunny ゲート付き往復テストで `energy(sys) ≈ predict_energy - j0` が
      等方・異方（DM/Γ）両方の fixture で成立（≥20 ランダム単位スピン配置、
      `atol ≈ 1e-8` eV）。
- [ ] `magesty sunny script <fixture>.xml -o out.jl` が `Meta.parseall` を通る
      スクリプトを出力し、`fept_tetragonal_2x2x2` で `dispersion` がエラーなく走り、
      Γ–X 経路に正の振動数ブランチが得られる。
- [ ] `make test-jet` / `make test-aqua` がクリーン。
- [ ] `docs/src/sunny_export.md`（使い方＋変換表）追加、`docs/make.jl` ナビ更新、
      CLI ドキュメントに追記。

## References

- 変換式: `docs/src/technical_notes.md`
- 既存パターン: `src/VaspConvert.jl`, `cli/src/MagestyCLI.jl`
- 関連 spec: 260520-cli-foundation, 260521-cli-package-extraction
- Sunny.jl: https://sunnysuite.github.io/Sunny.jl/stable/
