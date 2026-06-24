# Requirements: sce_to_sunny に物理スピン（実効スピン長）を渡す

Status: draft (2026-06-18)

## Goal

`sce_to_sunny` に**実効スピン長 `S_eff`**（= 角運動量 `ħS_eff`、`S_eff = m/(gμ_B)`）を
渡せるようにし、生成される Sunny スクリプトが物理的に正しいマグノン分散を出すよう
にする。現状は `Moment(s = 1)` 固定のため、SCE 係数が `J^SCE = J_phys·S²` の形で
スピンの大きさを吸収していることと相まって、LSWT マグノン振動数が真値の `S_eff`
倍だけ大きく出る（MnTe: Mn²⁺ S=5/2 で約 2.5 倍）。

## Background

- SCE は単位ベクトル `ê`（`‖ê‖=1`）でエネルギーをフィットするため、係数は
  `M^SCE = J_phys·S²` の形でモーメントの大きさを吸収している。静的エネルギーの
  再現には十分だが、**動力学（マグノン歳差）は角運動量 `ħS` に依存する**。
- LSWT マグノンは `ħω ∝ s·J_Sunny`。古典エネルギーを保つ再スケール
  `J_Sunny = M^SCE/s²` の下で `ħω ∝ M^SCE/s`。よって `s=1` 固定は真の `S_eff` に
  対して振動数を `S_eff` 倍に膨らませる。
- 実測値（INS）と DFT→SCE の整合は「同じ `S_eff` を両方向で使えば `J_phys` レベル
  で一致する」という論理。MnTe で数値確認済み: `M^SCE(最近接) = 24.5 meV`,
  `24.5/6.25 = 3.9 meV ≈ J₁ = 3.99 meV`（Liu et al., PRL 133, 156702 (2024)）。
- `S_eff` は物理的には半整数とは限らない（遍歴磁性 Fe ~2.2 μ_B → `S_eff ≈ 1.1`）が、
  **実装時に判明: Sunny の `Moment` は `s` を 1/2 の整数倍に強制**する
  （`"Spin s must be an exact multiple of 1/2"`、System 構築時に拒否）。よって
  非半整数 `S_eff` はそのまま渡せない。**決定（ユーザ承認, option B）: `spin` は
  半整数のみ許可し、非半整数は生成時に明確なエラーにする**。遍歴磁性対応
  （物理モーメントと Sunny の `s` を分離して分散だけ物理的にする option A）は将来の
  拡張として保留。
- 要求元: ユーザ（MnTe スピン波解析で 2 倍問題に直面）。

## Scope

Includes:

- `sce_to_sunny`（`src/SunnyExport.jl`）への公開キーワード追加:
  `spin`（必須, スカラ or 元素名→値の Dict）、`g`（任意, 既定 2）、
  `mode`（任意, `:auto`/`:dipole`/`:dipole_uncorrected`）。
- ボンドの再スケール `J_Sunny = M^SCE/(s_i·s_j)`、単一イオンの係数の
  `s_i`・モード依存化。
- 生成スクリプトの `Moment(s = …, g = …)` を副格子（元素）ごとの値に。
- CLI `magesty sunny script` に `--spin` ほかフラグを追加。
- ドキュメント（`docs/src/tools.md` の `magesty sunny script` 節 ＋ `api.md` の
  `@docs`）に `S_eff = m/(gμ_B)`、`s(2s-1)/2` リノーマライズ、モード選択を追記。

Excludes:

- SCE フィット側・XML 形式・SALC 構築の変更（一切触らない）。
- 外部磁場・中性子強度（フォームファクタ）の精密化（`g` は素通しのみ）。
- 縦ゆらぎ・Stoner 連続体など遍歴磁性の物理的限界そのものへの対処
  （ドキュメントで注意喚起するのみ）。

## Invariants

- **生成スクリプトの `energy(sys)` は `predict_energy(model) - j0` に一致**
  （任意の `spin`・`mode` で。`s_i s_j` が相殺するため古典エネルギーは spin 選択に
  不変）。これが正当性のアンカー。
- 分解 `M^SCE`（`_sunny_build_*` が返す素テンソル）と `predict_energy` の対応は
  不変。spin 再スケールは**emission 時にのみ**適用する。
- Spherical-harmonics 規約・SCE XML 形式・SALC↔Fitting の連結部は不変（読み取りのみ）。
- スピン方向は単位ベクトル・`3 × n_atoms` レイアウトのまま。

## Completion criteria

- [ ] `spin` 未指定時に物理的説明を含む `ArgumentError`（`UndefKeywordError` では
      なく、`S_eff = m/(gμ_B)` を案内する文言）。
- [ ] スカラ／Dict 両対応。非磁性種（couplings なし）は中立 placeholder で
      `Moment` を発行し、その旨をコメント。
- [ ] `energy(sys) ≈ predict_energy - j0` が `spin ∈ {1, 5/2, 1.1}` × モードで
      機械精度一致（Sunny ゲート付き往復テスト）。
- [ ] マグノン分散が `s=1` 基準に対し `1/s` でスケールする回帰テスト
      （バンド幅比 ≈ `s₁/s₂`）。
- [ ] `spin` が非半整数 → 生成時に `ArgumentError`（Sunny の制約を案内）。
- [ ] `mode = :auto` は半整数 `s` で `:dipole`。`:dipole_uncorrected` は明示指定で
      利用可（古典極限）。`:dipole` × `s=1/2` はエラー（四極子が消えるため）。
- [ ] `make test-all` / `make test-jet` / `make test-aqua` クリーン、かつ
      `make test-sunny`（Sunny ゲート、往復・分散スケール回帰）クリーン。
- [ ] `docs/src/tools.md`・`api.md`・CLI ヘルプ・`CHANGELOG.md`
      （`BREAKING CHANGE`）更新。

## References

- 関連 spec: [260529-sce-sunny-export](../260529-sce-sunny-export/)（本機能の土台）。
- 文献:
  - Z. Liu et al., "Chiral Split Magnon in Altermagnetic MnTe",
    Phys. Rev. Lett. **133**, 156702 (2024).（S=5/2、`H = Σ J S·S`、SpinW LSWT）
  - D. Dahlbom et al., "Renormalized Classical Theory of Quantum Magnets",
    arXiv:2304.03874.（`:dipole` の `c_k` リノーマライズ導出）
  - H. Zhang & C. D. Batista, "Classical spin dynamics based on SU(N) coherent
    states", Phys. Rev. B **104**, 104409 (2021).（古典スピン力学の基礎）
  - Sunny ドキュメント "Interaction Renormalization"
    （rank-2: `c₂ = 1 - 1/(2s)`）。
