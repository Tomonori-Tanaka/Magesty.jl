# Design: sce_to_sunny に物理スピン（実効スピン長）を渡す

Status: draft (2026-06-18)

## Summary

`sce_to_sunny` のシグネチャに `spin`（必須）・`g`・`mode` を追加する。分解段
（`_sunny_build_terms` / `_sunny_build_primitive`）は**現状のまま** `M^SCE` の素
テンソルを返し、**spin 再スケールは emission 段で適用**する。これにより:

- 分解と `predict_energy` の対応は不変（既存往復テストのオラクルを維持）。
- 生成スクリプトの `energy(sys) == predict_energy - j0` も任意の spin で不変
  （ボンドで `s_i s_j` が、単一イオンでモード補正が相殺するため）。
- 変更は emission 関数と新ヘルパ（spin 解決・onsite 係数）に局所化。

物理的背景（ドキュメント・コメントの根拠、ここには要約のみ）:

- マグノン歳差: `ħω ~ (角度に対するエネルギー曲率)/(角運動量 ħS)`。分子は
  `M^SCE`（DFT が固定、S 非依存）、分母は `ħS_eff = ħ·m/(gμ_B)`。よって
  `J_Sunny = M^SCE/(s_i s_j)`、`Moment(s = S_eff)` で物理マグノンを得る。
- `:dipole` の単一イオン補正 `s(2s-1)/2`: スピンコヒーレント状態の横ゆらぎ
  `⟨S_⊥²⟩ = s/2` が四極子を丸める量子効果。`s→∞` で `s²`（古典）、`s=1/2` で 0
  （四極子を持てない）。Sunny rank-2 係数 `c₂ = (2s-1)/(2s)`、`s(2s-1)/2 = c₂·s²`。

## Module layout

| Target | Change |
|---|---|
| `src/SunnyExport.jl` | `sce_to_sunny` に `spin`/`g`/`mode` 追加。`_SUNNY_SPIN_S` 定数を撤去し、emission 時に副格子ごとの `s` を解決して適用。新ヘルパ `_sunny_resolve_spins` / `_sunny_onsite_factor(s, mode)` / `_sunny_select_mode`。emit 関数（primitive/explicit）でボンドを `/(s_i s_j)`、onsite を係数化、`Moment` 行を per-site 化。ヘッダコメント更新。 |
| `cli/src/MagestyCLI.jl` | `sunny script` に `--spin`（`2.5` or `Mn=2.5,Fe=1.1`）・`--g`・`--mode` を追加し `sce_to_sunny` へ転送。 |
| `test/component/test_sunny_export.jl` | spin 解決・再スケール・onsite 係数・モード自動選択・エラー系の単体テスト追加（既存の `@test_throws ArgumentError ... placement=:nonsense` は**残したまま** `spin` 未指定のエラー系を追加）。 |
| `test/sunny/runtests.jl`（Sunny ゲート, `make test-sunny`。`test-all` には含まれない） | `spin × mode` の energy 不変性、`1/s` 分散スケール回帰を追加。新規ディレクトリは作らず既存ファイルに追記。 |
| `docs/src/tools.md`（`magesty sunny script` 節）, `docs/src/api.md`（`@docs sce_to_sunny`） | `S_eff` の意味・選び方、`s(2s-1)/2`、モード、遍歴磁性の注意、文献を追記。新規 `sunny_export.md` は作らない（`docs/make.jl` の `pages` 変更を避ける）。 |
| `CHANGELOG.md` | `[Unreleased]` に `BREAKING CHANGE`（`spin` 必須化）。 |

## API

```julia
function sce_to_sunny(
    model::SCEModel;
    spin::Union{Real, AbstractDict, Nothing} = nothing,  # REQUIRED. S_eff = m/(gμ_B)
    g::Union{Real, AbstractDict} = 2,
    mode::Symbol = :auto,                # :auto | :dipole | :dipole_uncorrected
    output::Union{AbstractString, Nothing} = nothing,
    placement::Symbol = :auto,           # 既存どおり
)::String
```

- `spin`: スカラ（全磁性種共通）または元素名 `String → Real` の Dict。物理スピン長
  `S_eff = m/(gμ_B)`。半整数に限らない。**既定 `nothing` は「有効な動作値」ではなく、
  `UndefKeywordError` ではなく `S_eff = m/(gμ_B)` を案内する `ArgumentError` を出す
  ためだけ**の番兵。実質必須。
- `g`: スカラ or 元素名 Dict。`Moment(s, g)` に素通し。**分散スケールには無影響**
  （外部磁場・強度向け）。
- `mode`: `:auto` は `:dipole`（spin は半整数に制約されるため）。`:dipole_uncorrected`
  （古典極限）は明示指定で利用可。`:dipole` は単一イオンに `s(2s-1)/2`、
  `:dipole_uncorrected` は `s²` を仮定。
- **Sunny 制約（実装時に確認）**: Sunny の `Moment` は `s` を 1/2 の整数倍に強制。
  非半整数 `spin` は `_sunny_resolve_spin_maps` で `ArgumentError`（壊れたスクリプト
  を出さない）。遍歴磁性は option B 採択により保留。

### 再スケールの規約（emission 段）

- **ボンド** `i–j`: `J_Sunny = M^SCE / (s_i · s_j)`（両モード共通、bilinear は
  サイト間でリノーマライズ無し）。
- **単一イオン** サイト `i`: 出力演算子の係数を
  - `:dipole` → `2/(s_i(2s_i-1))`（既存式の per-site 化。`= 1/(c₂ s²)`）
  - `:dipole_uncorrected` → `1/s_i²`（古典極限）
  - `s_i = 1/2` かつ `:dipole` かつ onsite 項あり → エラー（四極子が消え係数が発散）。
- **Moment 行**: 副格子（元素）ごとに解決した `s`・`g` を出力。
- **非磁性種**（couplings を持たない元素, 例 Te）: `Moment` 発行のため中立値
  （`s = 1`）を placeholder で割当て、磁気的に不活性である旨をコメント。分散に無影響。

### spin/g の解決（`_sunny_resolve_spins`）

- 入力がスカラ → 全副格子に同値。Dict → 元素名で照合、未掲載の元素は非磁性
  placeholder（`s=1`）。
- primitive ルート: 副格子の元素は `pm.prim.types`。explicit ルート: スーパーセル
  原子 `a` の元素は `struc.kd_name[struc.supercell.kd_int_list[a]]`（`kd_name[a]` の
  直接索引ではない。primitive ルート line 295 と同じ索引チェーン）。両ルートで同じ
  解決ロジックを共有。
- 磁性種（l≥1 の basis を持つ＝ couplings に登場する元素）が `spin` Dict に
  欠けていたらエラー（取りこぼし防止）。

## Types and conventions

- 数値規約変更: **`s=1` 既定の撤廃（必須化）は意図的な破壊的変更**。生成スクリプト
  のマグノン分散が変わる（静的エネルギーは不変）。理由・回帰テストを伴う。
- `energy(sys) == predict_energy - j0` の不変性は spin 選択・モードに依らず保持
  （`s_i s_j` と onsite 補正が相殺）。これを往復テストで `spin × mode` 横断で固定。
- onsite 係数の物理（`s(2s-1)/2`）はコード内コメントで自己完結に説明（外部
  scaffolding を参照しない）。
- 既存の引数なし `_sunny_onsite_factor()`（global `_SUNNY_SPIN_S` 依存, 定義 line 49 /
  使用 line 423）は `(s, mode)` 版で**完全に置換**。旧引数なし形は残さない。
- `_sunny_reconstruct_energy` と、それを使う Sunny 非依存テストは**未スケールの
  `M^SCE`** を検証するため**更新しない**。spin 再スケールは emission 段のみで、
  分解側（`_sunny_build_*`）には絶対に入れない（実装者の取り違え防止）。

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): 影響なし（変換は既存の
      `M0`/`A0` を再利用、読み取りのみ）。
- ~~SCE coefficient XML (`save` / `load`)~~: 影響なし。
- ~~`Fitting` ↔ `SALCBasis`~~: 影響なし（`_cluster_scaling` 等は読み取りのみ）。
- [ ] `SPEC.md` / `docs/src/api.md`: `sce_to_sunny` シグネチャ更新を反映。
      narrative ページ `docs/src/tools.md`（`magesty sunny script` 節）も併せて更新
      （api.md だけにしない）。

（`.claude/agents/` の sweep はモジュール名・Makefile 不変のため不要。exit
checklist に既出のためここでは除外。）

## Test strategy

- **単体（Sunny 非依存）** `test/component/test_sunny_export.jl`:
  - `_sunny_resolve_spins`: スカラ→全副格子、Dict→元素照合、非磁性 placeholder、
    磁性種欠落エラー、`spin=nothing` の `ArgumentError`。
  - 再スケール: ボンド行列が `M^SCE/(s_i s_j)` になる（リテラル値照合）。
  - onsite 係数: `:dipole`→`2/(s(2s-1))`、`:dipole_uncorrected`→`1/s²`、
    `s=1/2`+onsite+`:dipole`→エラー。
  - `_sunny_select_mode`: 半整数集合→`:dipole`、非半整数含む→`:dipole_uncorrected`。
  - `Moment(s, g)` 行の per-site 値。文字列の `Meta.parseall` 構文妥当性。
- **往復（Sunny ゲート付き）** `test/sunny/runtests.jl`（`make test-sunny`, `test-all`
  には非包含）:
  - `energy(sys) ≈ predict_energy - j0` を `spin ∈ {1, 5/2, 1.1}` × `mode` で機械
    精度確認（既存 fixture 流用、onsite ありの fixture を含める）。
  - **`1/s` 分散スケール**: 同一 fixture で `s=1` と `s=5/2` のバンド幅比 ≈ 2.5。
- **bench**: SunnyExport は hot path ではない（`Fitting`/`SALCBases` が hot）。
  `bench_log.md` 追記不要。

## Risks and open items

1. **破壊的変更**（`spin` 必須化）: 既存呼び出し・CLI・docs・examples が壊れる。
   全呼び出し箇所を洗い出し、`BREAKING CHANGE` を明記。
2. **`s=1/2` × onsite × `:dipole`**: 物理的に四極子が無い／係数発散。明示エラーで
   ガードし、`:dipole_uncorrected` を案内。
3. **非磁性種の Moment placeholder**: 不活性サイトが LSWT にゼロモードとして現れる
   可能性。既存 `s=1` 挙動を踏襲（現行スクリプトで実績あり）、回帰で確認。
4. **Sunny の半整数制約 → 遍歴磁性は保留（option B）**: Sunny は `s` を 1/2 の
   整数倍に強制するため、非半整数 `S_eff`（Fe 1.1 等）はそのまま渡せず、生成時に
   エラーとする。物理モーメントと Sunny の `s` を分離して分散だけ物理的にする
   option A（`J = M/(s_ref·S_eff)`、静的エネルギー不一致・単一イオン非対応）は
   将来拡張として保留。なお遍歴磁性自体の物理的限界（縦ゆらぎ・Stoner 連続体）も
   docs で注意喚起。
5. **`g` の役割**: 分散には無影響と明記（誤期待の防止）。

### 実装順（tasklist と対応）

spec 合意 → spin/g 解決ヘルパ＋モード選択 → emission 再スケール（primitive→explicit）
→ ヘッダ／Moment 行 → 単体テスト → Sunny 往復＋分散スケール回帰 → CLI → docs →
CHANGELOG。
